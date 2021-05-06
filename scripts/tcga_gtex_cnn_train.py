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
from keras.preprocessing.image import ImageDataGenerator
from keras.regularizers import l2

# trying to prevent cuDNN errors. taken from
# https://forums.developer.nvidia.com/t/could-not-create-cudnn-handle-cudnn-status-alloc-failed/108261/2
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

data = np.load(RESULT + "tcga_gtex_genes_data.npz")
x_train = data["x_train"]
x_test = data["x_test"]
y_train = data["y_train"].reshape(-1,1)
y_test = data["y_test"].reshape(-1,1)
y_train = keras.utils.to_categorical(y_train, 2)
y_test = keras.utils.to_categorical(y_test, 2)

print("x_train shape: ",   x_train.shape)
print("x_test shape:  ",    x_test.shape)
print("y_train shape: ",   y_train.shape)
print("y_test shape:  ",    y_test.shape)

model = Sequential(name="Sequential_" + dt_string)

model.add(Conv2D(96, (3, 3), activation='relu', padding = 'same', input_shape=(32, 32, 3)))    
model.add(Dropout(0.2))

model.add(Conv2D(96, (3, 3), activation='relu', padding = 'same'))  
model.add(Conv2D(96, (3, 3), activation='relu', padding = 'same', strides = 2))    
model.add(Dropout(0.5))

model.add(Conv2D(192, (3, 3), activation='relu', padding = 'same'))    
model.add(Conv2D(192, (3, 3), activation='relu', padding = 'same'))
model.add(Conv2D(192, (3, 3), activation='relu', padding = 'same', strides = 2))    
model.add(Dropout(0.5))    

model.add(Conv2D(192, (3, 3), padding = 'same'))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Conv2D(192, (1, 1),padding='valid'))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Conv2D(2, (1, 1), padding='valid'))

model.add(GlobalAveragePooling2D())

model.add(Activation('sigmoid'))
opt = optimizers.Adamax(lr=.0001)
model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
model.summary()

# Open the file
with open(MODELFOLDER + 'model_summary.txt','w') as fh:
    # Pass the file handle in as a lambda function to make it callable
    model.summary(print_fn=lambda x: fh.write(x + '\n'))

# Save the best model during each training checkpoint
checkpoint = ModelCheckpoint('.',
                            monitor='val_loss', 
                            verbose=0,
                            save_best_only= True,
                            mode='auto')

history = model.fit(x_train, y_train,callbacks=checkpoint,
                    epochs=40,
                    verbose=1,
                    validation_data = (x_test, y_test))

model.save(MODELFOLDER + "tcgamodel.h5", save_format='h5')

# summarize history for accuracy
plt.plot(history.history['accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.savefig(FIGURES + "accuracy.svg", dpi=300)
plt.close()

# summarize history for loss
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.savefig(FIGURES + "loss.svg", dpi=300)
plt.close()

# x_test has 3577 samples, up to 1500 0, afterwards 1
np.random.seed(1)
normal_idx = np.random.randint(1, 1500, size=4)
tumor_idx = normal_idx + 2000

#random_idx=range(4)
fig, axes = plt.subplots(1, 4, figsize=(16, 12))

for idx, ax in enumerate(axes.ravel()):
    img = x_test[normal_idx[idx]]
    #ax.set_title(y_test[idx])
    ax.axis('off')
    ax.imshow(img, aspect='equal')

#plt.gca().set_axis_off()
plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0.1)
plt.margins(5,5)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
#fig.suptitle('Normal Samples', fontsize=16)
plt.savefig(FIGURES + "normal.png", bbox_inches = 'tight', pad_inches = 0, dpi=300)

fig, axes = plt.subplots(1, 4, figsize=(16, 12))

for idx, ax in enumerate(axes.ravel()):
    img = x_test[tumor_idx[idx]]
    #ax.set_title(y_test[idx])
    ax.axis('off')
    ax.imshow(img, aspect='equal')

plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0.1)
plt.margins(5,5)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
#fig.suptitle('Normal Samples', fontsize=16)
plt.savefig(FIGURES + "tumor.png", bbox_inches = 'tight', pad_inches = 0, dpi=300)
