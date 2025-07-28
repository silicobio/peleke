from transformers import AutoTokenizer, AutoModelForCausalLM
from peft import PeftModel
import torch
import argparse
import warnings
import gc

# Suppress warnings
warnings.filterwarnings("ignore")

def generate_antibody_sequence(
    antigen_sequence: str,
    base_model_path: str = "microsoft/phi-4",
    adapter_path: str = "/home/nicholas/Documents/GitHub/peleke/models/peleke-phi-4/checkpoint-1914",
    max_new_tokens: int = 200,
    top_p: float = 0.9,
    temperature: float = 0.9,
    top_k: int = 50
) -> str:
    
    print("Loading tokenizer...")
    tokenizer = AutoTokenizer.from_pretrained(adapter_path, trust_remote_code=True)
    
    print("Loading base model...")
    # Load base model to CPU first to avoid device mapping issues
    base_model = AutoModelForCausalLM.from_pretrained(
        base_model_path,
        torch_dtype=torch.float16,
        trust_remote_code=True,
        device_map="cpu"  # Load to CPU first
    )
    
    print("Loading PEFT adapter...")
    model = PeftModel.from_pretrained(base_model, adapter_path)
    
    print("Moving to GPU...")
    # Move to GPU manually
    model = model.cuda()
    model.eval()
    
    # Direct generation without pipeline
    prompt = f"Antigen: {antigen_sequence}\nAntibody:"
    inputs = tokenizer(prompt, return_tensors="pt", truncation=True, max_length=400)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    
    print("Generating antibody...")
    with torch.no_grad():
        outputs = model.generate(
            **inputs,
            max_new_tokens=max_new_tokens,
            do_sample=True,
            top_p=top_p,
            temperature=temperature,
            top_k=top_k,
            pad_token_id=tokenizer.eos_token_id or tokenizer.pad_token_id,
            eos_token_id=tokenizer.eos_token_id
        )
    
    generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    
    # Clean up memory
    del model, base_model, inputs, outputs
    torch.cuda.empty_cache()
    gc.collect()
    
    return generated_text

def generate_antibody_sequence_conservative(
    antigen_sequence: str,
    base_model_path: str = "microsoft/phi-4",
    adapter_path: str = "/home/nicholas/Documents/GitHub/peleke/models/peleke-phi-4/checkpoint-1914",
    max_new_tokens: int = 200,
    temperature: float = 0.8
) -> str:
    """More conservative version with better memory management"""
    
    print("🔄 Conservative loading...")
    
    # Clear memory first
    torch.cuda.empty_cache()
    gc.collect()
    
    # Truncate very long antigens
    if len(antigen_sequence) > 350:
        print(f"⚠️  Truncating antigen from {len(antigen_sequence)} to 350 characters")
        antigen_sequence = antigen_sequence[:350]
    
    tokenizer = AutoTokenizer.from_pretrained(adapter_path, trust_remote_code=True)
    
    # Load with more conservative settings
    base_model = AutoModelForCausalLM.from_pretrained(
        base_model_path,
        torch_dtype=torch.float16,
        trust_remote_code=True,
        low_cpu_mem_usage=True
    )
    
    model = PeftModel.from_pretrained(base_model, adapter_path)
    model = model.cuda()
    model.eval()
    
    # Generate
    prompt = f"Antigen: {antigen_sequence}\nAntibody:"
    inputs = tokenizer(prompt, return_tensors="pt", truncation=True, max_length=350)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    
    with torch.no_grad():
        outputs = model.generate(
            **inputs,
            max_new_tokens=max_new_tokens,
            temperature=temperature,
            do_sample=True,
            top_p=0.9,
            pad_token_id=tokenizer.eos_token_id or tokenizer.pad_token_id,
            use_cache=True
        )
    
    generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    
    # Cleanup
    del model, base_model, inputs, outputs
    torch.cuda.empty_cache()
    gc.collect()
    
    return generated_text

def main():
    parser = argparse.ArgumentParser(description="Peleke🦋: Generate antibody sequences from antigen sequences.")
    parser.add_argument("--antigen", type=str, required=True,
                       help="The input antigen sequence with epitope residues surrounded by [square brackets].")
    parser.add_argument("--base_model_path", type=str, default="microsoft/phi-4",
                       help="Path to the base model (Phi-4).")
    parser.add_argument("--adapter_path", type=str,
                       default="/home/nicholas/Documents/GitHub/peleke/models/peleke-phi-4/checkpoint-1914",
                       help="Path to the fine-tuned PEFT adapter.")
    parser.add_argument("--max_new_tokens", type=int, default=200,
                       help="Maximum number of new tokens to generate.")
    parser.add_argument("--top_p", type=float, default=0.9,
                       help="Top-p sampling parameter.")
    parser.add_argument("--temperature", type=float, default=0.9,
                       help="Temperature for sampling.")
    parser.add_argument("--top_k", type=int, default=50,
                       help="Top-k sampling parameter.")
    parser.add_argument("--method", type=str, default="conservative", 
                       choices=["normal", "conservative"],
                       help="Generation method")
    
    args = parser.parse_args()
    
    print("🦋 Peleke Antibody Generator")
    print("=" * 50)
    print(f"🧬 Antigen length: {len(args.antigen)} characters")
    print(f"🔧 Method: {args.method}")
    
    try:
        if args.method == "conservative":
            antibody_sequence = generate_antibody_sequence_conservative(
                args.antigen,
                base_model_path=args.base_model_path,
                adapter_path=args.adapter_path,
                max_new_tokens=args.max_new_tokens,
                temperature=args.temperature
            )
        else:
            antibody_sequence = generate_antibody_sequence(
                args.antigen,
                base_model_path=args.base_model_path,
                adapter_path=args.adapter_path,
                max_new_tokens=args.max_new_tokens,
                top_p=args.top_p,
                temperature=args.temperature,
                top_k=args.top_k
            )
        
        print("\n" + "="*80)
        print("🎯 RESULT:")
        print("="*80)
        print(antibody_sequence)
        print("="*80)
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        
        if "compilation" in str(e).lower() or "python.h" in str(e).lower():
            print("\n💡 Compilation error detected. Install Python headers:")
            print("   sudo apt-get install python3.11-dev python3-dev build-essential")
        elif "memory" in str(e).lower() or "killed" in str(e).lower():
            print("\n💡 Memory error detected. Try:")
            print("   python script.py --method conservative --max_new_tokens 150")
        else:
            print(f"\n💡 Try the conservative method:")
            print("   python script.py --method conservative")

if __name__ == "__main__":
    main()