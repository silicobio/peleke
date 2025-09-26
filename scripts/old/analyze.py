import torch
import matplotlib.pyplot as plt
from transformers import AutoTokenizer
from peft import PeftModel

## Load model and tokenizer
model = PeftModel.from_pretrained("silicobio/peleke-phi4", device_map="auto", torch_dtype="auto")
tokenizer = AutoTokenizer.from_pretrained("silicobio/peleke-phi4", trust_remote_code=True)

## Example prompt
prompt = "Antigen: MKT[LLI]LAV[AA]A\nAntibody: "

inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
outputs = model(**inputs, output_attentions=True)
attentions = outputs.attentions[-1][0]  # Last layer, single head

# Get token positions
tokens = tokenizer.convert_ids_to_tokens(inputs["input_ids"][0])
epitope_indices = [i for i, t in enumerate(tokens) if "[" in t or "]" in t]
antibody_indices = [i for i, t in enumerate(tokens) if "|" in t or t.startswith("Q")]

# Plot attention from epitope tokens to antibody tokens
fig, ax = plt.subplots()
attn = attentions.mean(dim=0)[epitope_indices][:, antibody_indices]  # Mean over heads
cax = ax.matshow(attn.cpu().detach(), cmap="viridis")
plt.colorbar(cax)
plt.xticks(range(len(antibody_indices)), [tokens[i] for i in antibody_indices], rotation=90)
plt.yticks(range(len(epitope_indices)), [tokens[i] for i in epitope_indices])
plt.title("Epitope → Antibody Attention")
plt.tight_layout()
plt.show()
