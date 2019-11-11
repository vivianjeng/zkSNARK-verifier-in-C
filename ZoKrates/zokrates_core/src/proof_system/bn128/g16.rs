use crate::ir;
use crate::proof_system::bn128::utils::bellman::Computation;
use crate::proof_system::bn128::utils::solidity::{
    SOLIDITY_G2_ADDITION_LIB, SOLIDITY_PAIRING_LIB, SOLIDITY_PAIRING_LIB_V2,
};
use crate::proof_system::ProofSystem;
use bellman::groth16::Parameters;
use regex::Regex;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use zokrates_field::field::FieldPrime;

const G16_WARNING: &str = "WARNING: You are using the G16 scheme which is subject to malleability. See zokrates.github.io/reference/proving_schemes.html#g16-malleability for implications.";

pub struct G16 {}
impl ProofSystem for G16 {
    fn setup(&self, program: ir::Prog<FieldPrime>, pk_path: &str, vk_path: &str) {
        std::env::set_var("BELLMAN_VERBOSE", "0");

        println!("{}", G16_WARNING);

        let parameters = Computation::without_witness(program).setup();
        let parameters_file = File::create(PathBuf::from(pk_path)).unwrap();
        parameters.write(parameters_file).unwrap();
        let mut vk_file = File::create(PathBuf::from(vk_path)).unwrap();
        vk_file
            .write(serialize::serialize_vk(parameters.vk).as_ref())
            .unwrap();
    }

    fn generate_proof(
        &self,
        program: ir::Prog<FieldPrime>,
        witness: ir::Witness<FieldPrime>,
        pk_path: &str,
        proof_path: &str,
    ) -> bool {
        std::env::set_var("BELLMAN_VERBOSE", "0");

        println!("{}", G16_WARNING);

        let computation = Computation::with_witness(program, witness);
        let parameters_file = File::open(PathBuf::from(pk_path)).unwrap();

        let params = Parameters::read(parameters_file, true).unwrap();

        let proof = computation.clone().prove(&params);

        let mut proof_file = File::create(PathBuf::from(proof_path)).unwrap();
        write!(
            proof_file,
            "{}",
            serialize::serialize_proof(&proof, &computation.public_inputs_values())
        )
        .unwrap();
        true
    }

    fn export_solidity_verifier(&self, reader: BufReader<File>, is_abiv2: bool) -> String {
        let mut lines = reader.lines();

        let (mut template_text, solidity_pairing_lib) = if is_abiv2 {
            (
                String::from(CONTRACT_TEMPLATE_V2),
                String::from(SOLIDITY_PAIRING_LIB_V2),
            )
        } else {
            (
                String::from(CONTRACT_TEMPLATE),
                String::from(SOLIDITY_PAIRING_LIB),
            )
        };

        let gamma_abc_template = String::from("curvepoint_fp_init_set_str(vk_gamma_abc[index],points);"); //copy this for each entry

        //replace things in template
        let vk_regex = Regex::new(r#"(<%vk_[^i%]*%>)"#).unwrap();
        let vk_gamma_abc_len_regex = Regex::new(r#"(<%vk_gamma_abc_length%>)"#).unwrap();
        let vk_gamma_abc_index_regex = Regex::new(r#"index"#).unwrap();
        let vk_gamma_abc_points_regex = Regex::new(r#"points"#).unwrap();
        let vk_gamma_abc_repeat_regex = Regex::new(r#"(<%vk_gamma_abc_pts%>)"#).unwrap();
        let vk_input_len_regex = Regex::new(r#"(<%vk_input_length%>)"#).unwrap();

        for _ in 0..4 {
            let current_line: String = lines
                .next()
                .expect("Unexpected end of file in verification key!")
                .unwrap();
            let current_line_split: Vec<&str> = current_line.split("=").collect();
            assert_eq!(current_line_split.len(), 2);
            template_text = vk_regex
                .replace(template_text.as_str(), current_line_split[1].trim())
                .into_owned();
        }

        let current_line: String = lines
            .next()
            .expect("Unexpected end of file in verification key!")
            .unwrap();
        let current_line_split: Vec<&str> = current_line.split("=").collect();
        assert_eq!(current_line_split.len(), 2);
        let gamma_abc_count: i32 = current_line_split[1].trim().parse().unwrap();

        template_text = vk_gamma_abc_len_regex
            .replace(
                template_text.as_str(),
                format!("{}", gamma_abc_count).as_str(),
            )
            .into_owned();
        template_text = vk_input_len_regex
            .replace(
                template_text.as_str(),
                format!("{}", gamma_abc_count - 1).as_str(),
            )
            .into_owned();

        let mut gamma_abc_repeat_text = String::new();
        for x in 0..gamma_abc_count {
            let mut curr_template = gamma_abc_template.clone();
            let current_line: String = lines
                .next()
                .expect("Unexpected end of file in verification key!")
                .unwrap();
            let current_line_split: Vec<&str> = current_line.split("=").collect();
            assert_eq!(current_line_split.len(), 2);
            curr_template = vk_gamma_abc_index_regex
                .replace(curr_template.as_str(), format!("{}", x).as_str())
                .into_owned();
            curr_template = vk_gamma_abc_points_regex
                .replace(curr_template.as_str(), current_line_split[1].trim())
                .into_owned();
            gamma_abc_repeat_text.push_str(curr_template.as_str());
            if x < gamma_abc_count - 1 {
                gamma_abc_repeat_text.push_str("\n    ");
            }
        }
        template_text = vk_gamma_abc_repeat_regex
            .replace(template_text.as_str(), gamma_abc_repeat_text.as_str())
            .into_owned();

        let re = Regex::new(r"(?P<v>0[xX][0-9a-fA-F]{64})").unwrap();
        template_text = re.replace_all(&template_text, "$v").to_string();

        format!(
            "{}{}",
            solidity_pairing_lib, template_text
        )
    }
}

mod serialize {

    use crate::proof_system::bn128::utils::bellman::{
        parse_fr_json, parse_g1_hex, parse_g1_json, parse_g2_hex, parse_g2_json,
    };
    use bellman::groth16::{Proof, VerifyingKey};
    use pairing::bn256::{Bn256, Fr};

    pub fn serialize_vk(vk: VerifyingKey<Bn256>) -> String {
        format!(
            "vk.alpha = {}
    vk.beta = {}
    vk.gamma = {}
    vk.delta = {}
    vk.gamma_abc.len() = {}
    {}",
            parse_g1_hex(&vk.alpha_g1),
            parse_g2_hex(&vk.beta_g2),
            parse_g2_hex(&vk.gamma_g2),
            parse_g2_hex(&vk.delta_g2),
            vk.ic.len(),
            vk.ic
                .iter()
                .enumerate()
                .map(|(i, x)| format!("={}", parse_g1_hex(x)))
                .collect::<Vec<_>>()
                .join("\n")
        )
    }

    pub fn serialize_proof(p: &Proof<Bn256>, inputs: &Vec<Fr>) -> String {
        format!(
            "{} {} {} {}",
            parse_g1_json(&p.a),
            parse_g2_json(&p.b),
            parse_g1_json(&p.c),
            inputs
                .iter()
                .map(parse_fr_json)
                .collect::<Vec<_>>()
                .join(" "),
        )
    }
}

const CONTRACT_TEMPLATE_V2: &str = r#"
int main(int argc, char* argv[])
{
	// init
	fp_init();
	fp6_init();
	fp12_init();
	curve_init();
	points_init();

	// Input
	int num_elements = 8;
	int char_count = 79;
	int publicParamLeng = argc - 9;
	mpz_t publicParam[publicParamLeng];
	curvepoint_fp_t proof_A, proof_C;
	twistpoint_fp2_t proof_B;
	char **input = (char**) calloc(num_elements, sizeof(char*));
	for (int i=0;i< num_elements;i++)
	{
		input[i] = (char*) calloc(char_count, sizeof(char));
	}
	input[0] = argv[1];
	input[1] = argv[2];
	input[2] = argv[3];
	input[3] = argv[4];
	input[4] = argv[5];
	input[5] = argv[6];
	input[6] = argv[7];
	input[7] = argv[8];

	for(int i=0;i<publicParamLeng;i++){
		mpz_init_set_str(publicParam[i], argv[i+9], 10);
	}

	curvepoint_fp_init_set_str(proof_A, input[0],input[1],"1");
	twistpoint_fp2_init_set_str(proof_B,input[2],input[3],input[4],input[5]);
	curvepoint_fp_init_set_str(proof_C,input[6],input[7],"1");
	// verifying key
	curvepoint_fp_t vk_A;
	twistpoint_fp2_t vk_B, vk_gamma, vk_delta;
	curvepoint_fp_init_set_str(vk_A,<%vk_a%>); 
    twistpoint_fp2_init_set_str(vk_B,<%vk_b%>);
    twistpoint_fp2_init_set_str(vk_gamma,<%vk_gamma%>);
    twistpoint_fp2_init_set_str(vk_delta,<%vk_delta%>);
    curvepoint_fp_t vk_gamma_abc[<%vk_gamma_abc_length%>];
    <%vk_gamma_abc_pts%>

	// verify()
	curvepoint_fp_t vk_x, tmp0, tmp1, tmp2;
	curvepoint_fp_init_set_str(vk_x,"0","0","0");
	for(int i = 0; i< publicParamLeng; i++)
	{
		curvepoint_fp_mul(tmp0, vk_gamma_abc[i+1], publicParam[i]);
		curvepoint_fp_makeaffine(tmp0);
		if(fpe_iszero(vk_x->m_z)){
			curvepoint_fp_set(vk_x, tmp0);
		}
		else{
			curvepoint_fp_mixadd(vk_x, vk_x, tmp0);
			curvepoint_fp_makeaffine(vk_x);
		}
	}
	curvepoint_fp_mixadd(vk_x, vk_x, vk_gamma_abc[0]);
	curvepoint_fp_makeaffine(vk_x);

    // negate(vk_x)
    curvepoint_fp_neg(tmp0, vk_x);
	curvepoint_fp_makeaffine(tmp0);
    // negate(proof_C)
    curvepoint_fp_neg(tmp1, proof_C);
	curvepoint_fp_makeaffine(tmp1);
    // negate(vk_a)
    curvepoint_fp_neg(tmp2, vk_A);
	curvepoint_fp_makeaffine(tmp2);
	if(!pairingProd4(proof_A, proof_B, 
                     tmp0, vk_gamma,
                     tmp1, vk_delta,
                     tmp2, vk_B
                     )) return 1;
	
	printf("transaction verified");

    return 0;
}
"#;

const CONTRACT_TEMPLATE: &str = r#"
int main(int argc, char* argv[])
{
	// init
	fp_init();
	fp6_init();
	fp12_init();
	curve_init();
	points_init();

	// Input
	int num_elements = 8;
	int char_count = 79;
	int publicParamLeng = argc - 9;
	mpz_t publicParam[publicParamLeng];
	curvepoint_fp_t proof_A, proof_C;
	twistpoint_fp2_t proof_B;
	char **input = (char**) calloc(num_elements, sizeof(char*));
	for (int i=0;i< num_elements;i++)
	{
		input[i] = (char*) calloc(char_count, sizeof(char));
	}
	input[0] = argv[1];
	input[1] = argv[2];
	input[2] = argv[3];
	input[3] = argv[4];
	input[4] = argv[5];
	input[5] = argv[6];
	input[6] = argv[7];
	input[7] = argv[8];

	for(int i=0;i<publicParamLeng;i++){
		mpz_init_set_str(publicParam[i], argv[i+9], 10);
	}

	curvepoint_fp_init_set_str(proof_A, input[0],input[1],"1");
	twistpoint_fp2_init_set_str(proof_B,input[2],input[3],input[4],input[5]);
	curvepoint_fp_init_set_str(proof_C,input[6],input[7],"1");
	// verifying key
	curvepoint_fp_t vk_A;
	twistpoint_fp2_t vk_B, vk_gamma, vk_delta;
	curvepoint_fp_init_set_str(vk_A,<%vk_a%>); 
    twistpoint_fp2_init_set_str(vk_B,<%vk_b%>);
    twistpoint_fp2_init_set_str(vk_gamma,<%vk_gamma%>);
    twistpoint_fp2_init_set_str(vk_delta,<%vk_delta%>);
    curvepoint_fp_t vk_gamma_abc[<%vk_gamma_abc_length%>];
    <%vk_gamma_abc_pts%>

	// verify()
	curvepoint_fp_t vk_x, tmp0, tmp1, tmp2;
	curvepoint_fp_init_set_str(vk_x,"0","0","0");
	for(int i = 0; i< publicParamLeng; i++)
	{
		curvepoint_fp_mul(tmp0, vk_gamma_abc[i+1], publicParam[i]);
		curvepoint_fp_makeaffine(tmp0);
		if(fpe_iszero(vk_x->m_z)){
			curvepoint_fp_set(vk_x, tmp0);
		}
		else{
			curvepoint_fp_mixadd(vk_x, vk_x, tmp0);
			curvepoint_fp_makeaffine(vk_x);
		}
	}
	curvepoint_fp_mixadd(vk_x, vk_x, vk_gamma_abc[0]);
	curvepoint_fp_makeaffine(vk_x);

    // negate(vk_x)
    curvepoint_fp_neg(tmp0, vk_x);
	curvepoint_fp_makeaffine(tmp0);
    // negate(proof_C)
    curvepoint_fp_neg(tmp1, proof_C);
	curvepoint_fp_makeaffine(tmp1);
    // negate(vk_a)
    curvepoint_fp_neg(tmp2, vk_A);
	curvepoint_fp_makeaffine(tmp2);
	if(!pairingProd4(proof_A, proof_B, 
                     tmp0, vk_gamma,
                     tmp1, vk_delta,
                     tmp2, vk_B
                     )) return 1;
	
	printf("transaction verified");

    return 0;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;
    mod serialize {
        use super::*;

        mod proof {
            use super::*;
            use crate::flat_absy::FlatVariable;
            use crate::ir::*;
            use crate::proof_system::bn128::g16::serialize::serialize_proof;
            use typed_absy::types::{Signature, Type};

            #[allow(dead_code)]
            #[derive(Deserialize)]
            struct G16ProofPoints {
                a: [String; 2],
                b: [[String; 2]; 2],
                c: [String; 2],
            }

            #[allow(dead_code)]
            #[derive(Deserialize)]
            struct G16Proof {
                proof: G16ProofPoints,
                inputs: Vec<String>,
            }

            #[test]
            fn serialize() {
                let program: Prog<FieldPrime> = Prog {
                    main: Function {
                        id: String::from("main"),
                        arguments: vec![FlatVariable::new(0)],
                        returns: vec![FlatVariable::public(0)],
                        statements: vec![Statement::Constraint(
                            FlatVariable::new(0).into(),
                            FlatVariable::public(0).into(),
                        )],
                    },
                    private: vec![false],
                    signature: Signature::new()
                        .inputs(vec![Type::FieldElement])
                        .outputs(vec![Type::FieldElement]),
                };

                let witness = program
                    .clone()
                    .execute(&vec![FieldPrime::from(42)])
                    .unwrap();
                let computation = Computation::with_witness(program, witness);

                let public_inputs_values = computation.public_inputs_values();

                let params = computation.clone().setup();
                let proof = computation.prove(&params);

                let serialized_proof = serialize_proof(&proof, &public_inputs_values);
                serde_json::from_str::<G16Proof>(&serialized_proof).unwrap();
            }
        }
    }
}
