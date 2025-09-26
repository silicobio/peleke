from transformers import AutoTokenizer, AutoModelForCausalLM, pipeline
import torch
import argparse

def generate_antibody_sequence(
        antigen_sequence: str,
        model_path: str="../models/peleke-phi-4",
        max_new_tokens: int=200,
        top_p: float=0.9,
        temperature: float=0.9,
        top_k: int=50
        ) -> str:
    tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True)
    model = AutoModelForCausalLM.from_pretrained(model_path, device_map="auto", torch_dtype=torch.float16)

    generator = pipeline("text-generation", model=model, tokenizer=tokenizer, device=0)
    prompt = f"Antigen: {antigen_sequence}\nAntibody:"

    outputs = generator(prompt,
                        max_new_tokens=max_new_tokens,
                        do_sample=True,
                        top_p=top_p,
                        temperature=temperature,
                        top_k=top_k)

    return outputs[0]["generated_text"]


def main():
    parser = argparse.ArgumentParser(description="Peleke🦋: Generate antibody sequences from antigen sequences.")
    parser.add_argument("--antigen", type=str, required=True, help="The input antigen sequence with epitope residues surrounded by [square brackets].")
    parser.add_argument("--model_path", type=str, default="../models/peleke-phi-4", help="Path to the pre-trained model.")
    parser.add_argument("--max_new_tokens", type=int, default=200, help="Maximum number of new tokens to generate.")
    parser.add_argument("--top_p", type=float, default=0.9, help="Top-p sampling parameter.")
    parser.add_argument("--temperature", type=float, default=0.9, help="Temperature for sampling.")
    parser.add_argument("--top_k", type=int, default=50, help="Top-k sampling parameter.")
    
    args = parser.parse_args()

    antibody_sequence = generate_antibody_sequence(
        args.antigen,
        model_path=args.model_path,
        max_new_tokens=args.max_new_tokens,
        top_p=args.top_p,
        temperature=args.temperature,
        top_k=args.top_k
    )

    print(antibody_sequence)


if __name__ == "__main__":
    main()