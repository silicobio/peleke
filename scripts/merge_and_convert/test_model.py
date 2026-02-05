import torch, re
from transformers import AutoModelForCausalLM, AutoTokenizer

import os
os.environ["TRANSFORMERS_OFFLINE"] = "1"
os.environ["HF_HOME"] = "/scratch/cford38/peleke_tests/.cache/huggingface/"

# model_path = "./peleke-phi-4-merged"
model_path = "silicobio/peleke-phi-4-merged"

model = AutoModelForCausalLM.from_pretrained(
    model_path,
    torch_dtype=torch.bfloat16,
    device_map="auto",
    trust_remote_code=True,
)

tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True)

def format_prompt(antigen_sequence):
    epitope_seq = re.sub(r'\[([A-Z])\]', r'<epi>\1</epi>', antigen_sequence)
    formatted_str = f"Antigen: {epitope_seq}<|im_end|>\nAntibody:"
    return formatted_str

## Example HIV-1 gp120 antigen sequence with epitopes marked
antigen_sequence = "VPVWKDADTTLFCASDAKAHETEVHNVWATHACVPTDPNPQEIHLENVTENFNMWKNNMVEQMQEDVISLWDQSLQPCVKLTGGSVIKQACPKISFDPIPIHYCTPAGYVILKCNDKNFNGTGPCKNVSSVQCTHGIKPVVSTQLLLNGSLAEEEIIIRSENLT[N][N]A[K]TIIVHLNKSVEINCTRPSNGGSGSGGDIRKAYCEINGTKWNKVLKQVTEKLKEHFNNKTIIFQPPSGG[D]LEITMHSFNCRGEFFYCNTTQLFNNTCIGNETMKGCNGTITLPCKIKQIINMWQGTGQAMYAPPIDGKINCVSNITGILLT[R]D[G][G]ANNTSNETFRPGGGNIKDNWRSELYKYKVVQIE"

prompt = format_prompt(antigen_sequence)
inputs = tokenizer(prompt, return_tensors="pt")
inputs = {k: v.cuda() for k, v in inputs.items()}

with torch.no_grad():
    outputs = model.generate(
        **inputs,
        max_new_tokens=1000,
        do_sample=True,
        temperature=0.7,
        pad_token_id=tokenizer.eos_token_id,
        use_cache=False,
    )

full_text = tokenizer.decode(outputs[0], skip_special_tokens=False)
antibody_sequence = full_text.split('<|im_end|>')[1].replace('Antibody: ', '')
print(f"Antigen: {antigen_sequence.replace('[','').replace(']','')}\nAntibody: {antibody_sequence}\n")
