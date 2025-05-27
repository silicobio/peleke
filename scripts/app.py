import gradio as gr
from transformers import AutoTokenizer, AutoModelForCausalLM, pipeline

# Load model and tokenizer
model_name = "tuplexyz/peleke-phi4"
tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
model = AutoModelForCausalLM.from_pretrained(model_name, device_map="auto")

generator = pipeline("text-generation", model=model, tokenizer=tokenizer, device=0)

def generate_antibody(antigen):
    prompt = f"Antigen: {antigen}\nAntibody:"
    out = generator(prompt, max_new_tokens=200, temperature=0.9)[0]["generated_text"]
    return out.split("Antibody:")[-1].strip()

gr.Interface(
    fn=generate_antibody,
    inputs=gr.Textbox(label="Epitope-Annotated Antigen Sequence", placeholder="MKT[LLI]LAV[AA]A"),
    outputs=gr.Textbox(label="Generated Antibody Sequence"),
    title="Antibody Generator from Epitope",
    description="Enter an epitope-tagged antigen amino acid sequence to generate a heavy|light Fv antibody sequence."
).launch()
