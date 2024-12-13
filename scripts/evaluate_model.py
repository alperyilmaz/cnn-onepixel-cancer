# Python Libraries
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import backend
from tensorflow.keras.datasets import cifar10
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Conv2D, Dense, Flatten, MaxPooling2D, Dropout, Activation, GlobalAveragePooling2D
from tensorflow.keras.callbacks import LearningRateScheduler, TensorBoard, ModelCheckpoint
from tensorflow.keras.regularizers import l2
from PIL import Image
from datetime import datetime
from sklearn.metrics import roc_curve, auc, classification_report, confusion_matrix

# GPU Configuration
gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        logical_gpus = tf.config.list_logical_devices('GPU')
        print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
    except RuntimeError as e:
        print(e)

# Paths and Constants
RESULT = "results/"
MODELFOLDER = "model/"
FIGURES = "figures/"

now = datetime.now()
dt_string = now.strftime("%Y%m%d-%H%M%S")

class_names = ["Normal", "Tumor"]

# Load Data
data = np.load(RESULT + "tcga_gtex_genes_data.npz")
x_train = data["x_train"]
x_test = data["x_test"]
y_train = data["y_train"].reshape(-1, 1)
y_test = data["y_test"].reshape(-1, 1)

# Load Model and Predict
tcganet = keras.models.load_model(MODELFOLDER + "tcgamodel.h5")
predictions = tcganet.predict(x_test)
predictions_cold = [np.argmax(i) for i in predictions]

# Model Evaluation
with open(MODELFOLDER + "model_evaluation_report.txt", 'w') as f:
    print('Confusion Matrix', file=f)
    print(confusion_matrix(y_test[:, 0], predictions_cold), file=f)
    print("\nClassification Report", file=f)
    print(classification_report(y_test[:, 0], predictions_cold, target_names=class_names), file=f)

# ROC Curve
fpr_keras, tpr_keras, _ = roc_curve(y_test[:, 0], predictions[:, 1])
auc_keras = auc(fpr_keras, tpr_keras)

plt.figure(figsize=(8, 6))
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr_keras, tpr_keras, label=f'Area = {auc_keras:.3f}')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC Curve')
plt.legend(loc='best')
plt.savefig(FIGURES + "ROC_AUC.svg", dpi=300)
plt.close()
