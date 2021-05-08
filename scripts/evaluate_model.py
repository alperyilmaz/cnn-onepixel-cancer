# Python Libraries
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from keras.datasets import cifar10
from keras import backend as K
import tensorflow as tf
import keras
from PIL import Image
from datetime import datetime
from keras import optimizers
from keras.datasets import cifar10
from keras.models import Sequential, load_model
from keras.layers import Conv2D, Dense, Flatten, MaxPooling2D, Dropout, Activation, GlobalAveragePooling2D
from keras.callbacks import LearningRateScheduler, TensorBoard, ModelCheckpoint
from keras.regularizers import l2

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
  try:
    # Currently, memory growth needs to be the same across GPUs
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
    logical_gpus = tf.config.experimental.list_logical_devices('GPU')
    print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
  except RuntimeError as e:
    # Memory growth must be set before GPUs have been initialized
    print(e)

RESULT="results/"
MODELFOLDER="model/"
FIGURES="figures/"

now = datetime.now()
dt_string = now.strftime("%Y%m%d-%H%M%S")

class_names = ["Normal","Tumor"]
data = np.load(RESULT + "tcga_gtex_genes_data.npz")
x_train = data["x_train"]
x_test = data["x_test"]
y_train = data["y_train"].reshape(-1,1)
y_test = data["y_test"].reshape(-1,1)

## print more model related information
from sklearn.metrics import roc_curve, auc, roc_auc_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt2

tcganet = keras.models.load_model(MODELFOLDER + "tcgamodel.h5")
predictions = tcganet.predict(x_test)
true_pred_tuples = [(label,np.argmax(pred)) for (label, pred) in zip(y_test[:, 0], predictions)]
predictions_cold = [ np.argmax(i) for i in predictions ]

f=open(MODELFOLDER + "model_evaluation_report.txt",'w')
print('Confusion Matrix', file=f)
print(confusion_matrix(y_test[:, 0], predictions_cold), file=f)
print(" ",file=f)
print('Classification Report',file=f)
print(classification_report(y_test[:, 0], predictions_cold, target_names=class_names), file=f)
f.close()

fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test[:, 0], predictions_cold)
auc_keras = auc(fpr_keras, tpr_keras)

plt2.figure(1)
plt2.plot([0, 1], [0, 1], 'k--')
plt2.plot(fpr_keras, tpr_keras, label='area = {:.3f}'.format(auc_keras))
plt2.xlabel('False positive rate')
plt2.ylabel('True positive rate')
plt2.title('ROC curve')
plt2.legend(loc='best')
plt2.savefig(FIGURES + "ROC_AUC.svg", dpi=300)
plt2.close()
