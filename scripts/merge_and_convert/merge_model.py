import torch
from transformers import AutoModelForCausalLM, AutoTokenizer
from peft import PeftModel, PeftConfig

import os
os.environ["TRANSFORMERS_OFFLINE"] = "1"
os.environ["HF_HOME"] = "/scratch/cford38/peleke_tests/.cache/huggingface/"

## Parge Model name from Command Line
parser = argparse.ArgumentParser(description="Merging PEFT/LoRA models into singular models.")
parser.add_argument("--model_name", type=str, required=True, help="Model name from Hugging Face.")
args = parser.parse_args()

## Define model name and output directory
# model_name = "silicobio/peleke-phi-4"
model_name = args.model_name
output_dir = f"./{model_name.split('/')[1]}-merged"

## Load PEFT config to get base model
config = PeftConfig.from_pretrained(model_name)

## Load tokenizer (from PEFT repo is usually correct)
tokenizer = AutoTokenizer.from_pretrained(
    model_name,
    trust_remote_code=True
)

## Load base model
model = AutoModelForCausalLM.from_pretrained(
    config.base_model_name_or_path,
    torch_dtype=torch.bfloat16,
    trust_remote_code=True,
    device_map="auto",
)

## Make sure embeddings match tokenizer
model.resize_token_embeddings(len(tokenizer))

## Load PEFT weights on top
model = PeftModel.from_pretrained(
    model,
    model_name,
    torch_dtype=torch.bfloat16,
)

## Merge PEFT weights into the base model
model = model.merge_and_unload()

## Save merged model + tokenizer
model.save_pretrained(output_dir, safe_serialization=True)
tokenizer.save_pretrained(output_dir)
print(f"Merged model saved to: {output_dir}")


## Fix Tokenizer Mistral regex issue and resave
tokenizer = AutoTokenizer.from_pretrained(
    output_dir,
    use_fast=True,
    trust_remote_code=True,
    fix_mistral_regex=True,
)

tokenizer.save_pretrained(output_dir)
print(f"Merged tokenizer saved to: {output_dir}")