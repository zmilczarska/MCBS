# -*- coding: utf-8 -*-
"""pneumonia_classification_vit_homework.ipynb

Homework assignment: PneumoniaMNIST classification with Vision Transformer
"""
#pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128
#pip install numpy
#pip install matplotlib
#pip install medmnist
#pip install transformers

import sys
import os
import torch
from torch import nn
import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader
from torchvision import transforms
from medmnist import PneumoniaMNIST
from transformers import ViTConfig, ViTModel

# Random seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)

# Device configuration
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {DEVICE}")

# Dataset and model parameters
BATCH_SIZE = 32
LEARNING_RATE = 1e-4
NUM_EPOCHS = 10
IMAGE_SIZE = 224  # ViT typically uses 224x224 images
PATCH_SIZE = 16
NUM_CLASSES = 1  # Binary classification (normal vs pneumonia)

# Transforms
TRANSFORMS = transforms.Compose([
    transforms.ToTensor(),  # Converts PIL Image to tensor and scales to [0, 1]
    transforms.Lambda(lambda x: x.repeat(3, 1, 1)),  # Convert grayscale to 3 channels
    transforms.Resize((IMAGE_SIZE, IMAGE_SIZE)),
    transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])

# Load PneumoniaMNIST dataset
data_dir = "./pneumonia_data"
os.makedirs(data_dir, exist_ok=True)

train_set = PneumoniaMNIST(root=data_dir, split="train", transform=TRANSFORMS, download=True)
val_set = PneumoniaMNIST(root=data_dir, split="val", transform=TRANSFORMS)
test_set = PneumoniaMNIST(root=data_dir, split="test", transform=TRANSFORMS)

# DataLoaders
train_loader = DataLoader(train_set, batch_size=BATCH_SIZE, shuffle=True)
val_loader = DataLoader(val_set, batch_size=BATCH_SIZE, shuffle=False)
test_loader = DataLoader(test_set, batch_size=BATCH_SIZE, shuffle=False)

# Class names
CLASS_NAMES = ["normal", "pneumonia"]
LOGIT2NAME = {0: "normal", 1: "pneumonia"}

# Visualize some samples
samples, labels = next(iter(train_loader))
fig, axes = plt.subplots(1, 5, figsize=(15, 3))
for i in range(5):
    # Reverse normalization for visualization
    img = samples[i].permute(1, 2, 0).cpu().numpy()
    img = img * np.array([0.229, 0.224, 0.225]) + np.array([0.485, 0.456, 0.406])
    img = np.clip(img, 0, 1)
    axes[i].imshow(img)
    axes[i].set_title(f"Label: {CLASS_NAMES[labels[i].item()]}")
    axes[i].axis('off')
plt.tight_layout()
plt.savefig("samples.png")
plt.show()


# Define a simpler ViT model
class PneumoniaViT(nn.Module):
    def __init__(self):
        super().__init__()
        # Smaller ViT configuration than DINO
        config = ViTConfig(
            image_size=IMAGE_SIZE,
            patch_size=PATCH_SIZE,
            num_channels=3,
            hidden_size=192,  # Smaller than DINO's 384
            num_hidden_layers=6,  # Fewer layers than DINO's 12
            num_attention_heads=6,
            intermediate_size=768,
            num_labels=NUM_CLASSES,
            attention_probs_dropout_prob=0.1,
            hidden_dropout_prob=0.1
        )
        self.backbone = ViTModel(config)
        self.head = nn.Linear(config.hidden_size, NUM_CLASSES)

    def forward(self, x, output_attentions=False):
        outputs = self.backbone(x, output_attentions=output_attentions)
        x = outputs.last_hidden_state[:, 0]  # Use [CLS] token
        x = self.head(x)
        if output_attentions:
            return x, outputs.attentions
        return x


# Initialize model
model = PneumoniaViT().to(DEVICE)

# Loss and optimizer
criterion = nn.BCEWithLogitsLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE)

# Training loop
best_val_loss = float('inf')
for epoch in range(NUM_EPOCHS):
    model.train()
    train_loss = 0.0

    for images, labels in train_loader:
        images, labels = images.to(DEVICE), labels.float().view(-1, 1).to(DEVICE)

        # Forward pass
        outputs = model(images)
        loss = criterion(outputs, labels)

        # Backward pass and optimize
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        train_loss += loss.item()

    # Validation
    model.eval()
    val_loss = 0.0
    correct = 0
    total = 0

    with torch.no_grad():
        for images, labels in val_loader:
            images, labels = images.to(DEVICE), labels.float().view(-1, 1).to(DEVICE)
            outputs = model(images)
            loss = criterion(outputs, labels)
            val_loss += loss.item()

            preds = torch.sigmoid(outputs).round()
            correct += (preds == labels).sum().item()
            total += labels.size(0)

    train_loss /= len(train_loader)
    val_loss /= len(val_loader)
    val_acc = correct / total

    print(f"Epoch [{epoch + 1}/{NUM_EPOCHS}], "
          f"Train Loss: {train_loss:.4f}, "
          f"Val Loss: {val_loss:.4f}, "
          f"Val Acc: {val_acc:.4f}")

    # Save best model
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        torch.save(model.state_dict(), "best_pneumonia_vit.pth")

# Load best model
model.load_state_dict(torch.load("best_pneumonia_vit.pth", map_location=DEVICE))
model.eval()

# Test evaluation
correct = 0
total = 0
with torch.no_grad():
    for images, labels in test_loader:
        images, labels = images.to(DEVICE), labels.float().view(-1, 1).to(DEVICE)
        outputs = model(images)
        preds = torch.sigmoid(outputs).round()
        correct += (preds == labels).sum().item()
        total += labels.size(0)

print(f"Test Accuracy: {correct / total:.4f}")

# Explainability analysis
# Get some test samples
test_samples, test_labels = next(iter(test_loader))
test_samples = test_samples[:5].to(DEVICE)
test_labels = test_labels[:5].float().view(-1, 1).to(DEVICE)
print("test_samples shape:", test_samples.shape)
print("test_labels shape:", test_labels.shape)


# Predictions
with torch.no_grad():
    logits = model(test_samples)
    predictions = torch.sigmoid(logits).round().squeeze().cpu().numpy()

# Print predictions
for i, (pred, label) in enumerate(zip(predictions, test_labels)):
    print(f"Sample {i + 1}:")
    print(f"  True Label: {CLASS_NAMES[int(label.item())]}")
    print(f"  Prediction: {CLASS_NAMES[int(pred)]}")
    print(f"  Confidence: {torch.sigmoid(logits[i]).item():.4f}")

# Attention maps
with torch.no_grad():
    _, attentions = model(test_samples, output_attentions=True)

# Get attention from last layer and average heads
attention_maps = attentions[-1].mean(dim=1)[:, 0, 1:]
attention_maps = attention_maps.reshape(-1, IMAGE_SIZE // PATCH_SIZE, IMAGE_SIZE // PATCH_SIZE)

# Resize attention maps to match image size
attention_maps_resized = torch.nn.functional.interpolate(
    attention_maps.unsqueeze(1),
    size=(IMAGE_SIZE, IMAGE_SIZE),
    mode='bilinear',
    align_corners=False
).squeeze(1)

# Visualize attention
fig, axes = plt.subplots(2, 5, figsize=(15, 6))
for i in range(5):
    # Original image
    img = test_samples[i].permute(1, 2, 0).cpu().numpy()
    img = img * np.array([0.229, 0.224, 0.225]) + np.array([0.485, 0.456, 0.406])
    img = np.clip(img, 0, 1)
    axes[0, i].imshow(img)
    axes[0, i].set_title(f"{CLASS_NAMES[int(test_labels[i].item())]}")
    axes[0, i].axis('off')

    # Attention map
    attn = attention_maps_resized[i].cpu().numpy()
    axes[1, i].imshow(attn, cmap='hot')
    axes[1, i].set_title(f"Pred: {CLASS_NAMES[int(predictions[i])]}")
    axes[1, i].axis('off')

plt.tight_layout()
plt.savefig("attention_explanation.png")
plt.show()

# Gradient-based explanation
test_samples = test_samples.detach().requires_grad_(True)
logits = model(test_samples)
loss = criterion(logits, test_labels)
loss.backward()

# Get gradients and process them
grads = test_samples.grad.abs().mean(dim=1)  
grads = grads.detach().cpu().numpy()

# Visualize gradients
fig, axes = plt.subplots(1, 5, figsize=(15, 3))
for i in range(5):
    img = test_samples[i].permute(1, 2, 0).detach().cpu().numpy()
    img = img * np.array([0.229, 0.224, 0.225]) + np.array([0.485, 0.456, 0.406])
    img = np.clip(img, 0, 1)

    grad = grads[i]
    grad = (grad - grad.min()) / (grad.max() - grad.min())  # Normalize

    axes[i].imshow(img, cmap='gray')
    axes[i].imshow(grad, cmap='hot', alpha=0.5)
    axes[i].set_title(f"{CLASS_NAMES[int(predictions[i])]}")
    axes[i].axis('off')

plt.tight_layout()
plt.savefig("gradient_explanation.png")
plt.show()

print("All steps completed successfully.")
