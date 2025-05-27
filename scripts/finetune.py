from transformers import AutoTokenizer, AutoModelForCausalLM, Trainer, TrainingArguments
from datasets import load_dataset, Dataset
from peft import get_peft_model, LoraConfig, TaskType
import torch

import pandas as pd

## Load dataset
df = pd.read_csv("../data/train.csv")
df["text"] = df.apply(lambda x: f"Antigen: {x['antigen']}\nAntibody: {x['antibody']}", axis=1)
dataset = Dataset.from_pandas(df)

## Format prompts
def format_prompt(example):
    return {
        "text": f"Antigen: {example['antigen']}\nAntibody: {example['antibody']}"
    }

dataset = Dataset.from_pandas(df)
dataset = dataset.map(format_prompt)


## Load base tokenizer and model
model_name = "microsoft/phi-4"
tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
model = AutoModelForCausalLM.from_pretrained(model_name, device_map="auto", torch_dtype=torch.float16)

## Extend tokenizer with special tokens
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
extra_tokens = amino_acids + ["[", "]", "|"]

## Check if tokens already exist in the tokenizer's vocabulary
new_tokens = [t for t in extra_tokens if t not in tokenizer.get_vocab()]
tokenizer.add_tokens(new_tokens)
model.resize_token_embeddings(len(tokenizer))

## Tokenize the dataset
def tokenize(example):
    return tokenizer(example["text"], padding="max_length", truncation=True, max_length=512)

tokenized_dataset = dataset.map(tokenize)

## Training arguments
training_args = TrainingArguments(
    output_dir="./peleke-phi4",
    per_device_train_batch_size=2,
    per_device_eval_batch_size=2,
    num_train_epochs=3,
    warmup_steps=100,
    weight_decay=0.01,
    logging_dir="./logs",
    logging_steps=50,
    save_strategy="epoch",
    evaluation_strategy="no",
    fp16=True,
    # gradient_checkpointing=True ## If having memory issues
)

## PEFT configuration
peft_config = LoraConfig(
    r=8,
    lora_alpha=16,
    lora_dropout=0.05,
    bias="none",
    task_type=TaskType.CAUSAL_LM,
)

model = get_peft_model(model, peft_config)
model.print_trainable_parameters()

## Trainer
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_dataset,
    tokenizer=tokenizer,
)

## Fine-tune
trainer.train()
