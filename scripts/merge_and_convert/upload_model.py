from transformers import AutoModelForCausalLM, AutoTokenizer
import torch
import os

from huggingface_hub import login
login(token=os.environ["HF_TOKEN"])

model_name = "silicobio/peleke-phi-4-merged"
model_path = "./peleke-phi-4-merged"

model = AutoModelForCausalLM.from_pretrained(
    model_path,
    torch_dtype=torch.bfloat16,
    device_map="auto",
    trust_remote_code=True,
)
model.push_to_hub(model_name)

tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True)
tokenizer.push_to_hub(model_name)
