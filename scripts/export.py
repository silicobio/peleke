from transformers import AutoTokenizer
from peft import PeftModel, PeftConfig
from huggingface_hub import notebook_login

model = PeftModel.from_pretrained("tuplexyz/peleke-phi4", device_map="auto", torch_dtype="auto")
tokenizer = AutoTokenizer.from_pretrained("tuplexyz/peleke-phi4", trust_remote_code=True)

## Push LoRA adapters only
model.push_to_hub("tuplexyz/peleke-phi4", use_auth_token=True)
tokenizer.push_to_hub("tuplexyz/peleke-phi4", use_auth_token=True)
