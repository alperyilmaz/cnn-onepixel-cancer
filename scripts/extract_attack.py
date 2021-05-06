# Python Libraries
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from keras import backend as K
import tensorflow as tf
import keras

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
MODEL="model/"

tcganet = keras.models.load_model(MODEL + "tcgamodel.h5")
models = [tcganet]

class_names = ["Normal","Tumor"]
data = np.load(RESULT + "tcga_gtex_genes_data.npz", allow_pickle=True)
x_train = data["x_train"]
x_test = data["x_test"]
y_train = data["y_train"].reshape(-1,1)
y_test = data["y_test"].reshape(-1,1)
gene_names = data["genes"]
test_sample_names = data["test_samples"]
control= data["control"]

def predict_classes(xs, img, target_class, model, minimize=True):
    # Perturb the image with the given pixel(s) x and get the prediction of the model
    imgs_perturbed = perturb_image(xs, img)
    predictions = model.predict(imgs_perturbed)[:,target_class]
    # This function should always be minimized, so return its complement if needed
    return predictions if minimize else 1 - predictions

def color_process(imgs):
    if imgs.ndim < 4:
        imgs = np.array([imgs])
    imgs  = imgs.astype('float32')
    std   = [np.std(x_test[:,:,:,i]) for i in range(3)]
    mean  = [np.mean(x_test[:,:,:,i]) for i in range(3)]
    for img in imgs:
        for i in range(3):
            img[:,:,i] = (img[:,:,i] - mean[i]) / std[i]
    return imgs

def predict(img):
    processed = color_process(img)
    return tcganet.predict(processed, batch_size=64)

def predict_one(img):
    return predict(img)[0]

def attack_success(x, img, target_class, model, targeted_attack=False, verbose=False):
    # Perturb the image with the given pixel(s) and get the prediction of the model
    attack_image = perturb_image(x, img)

    confidence = model.predict(attack_image)[0]
    predicted_class = np.argmax(confidence)
    
    # If the prediction is what we want (misclassification or 
    # targeted classification), return True
    if verbose:
        print('Confidence:', confidence[target_class])
    if ((targeted_attack and predicted_class == target_class) or
        (not targeted_attack and predicted_class != target_class)):
        return True
    # NOTE: return None otherwise (not False), due to how Scipy handles its callback function

def attack(img_id, model, target=None, pixel_count=1, 
           maxiter=75, popsize=400, verbose=False):
    # Change the target class based on whether this is a targeted attack or not
    targeted_attack = target is not None
    target_class = target if targeted_attack else y_test[img_id, 0]
    
    # Define bounds for a flat vector of x,y,r,g,b values
    # For more pixels, repeat this layout
    bounds = [(0,32), (0,32), (0,2), (0,256), (0,256)] * pixel_count
    
    # Population multiplier, in terms of the size of the perturbation vector x
    popmul = max(1, popsize // len(bounds))
    
    # Format the predict/callback functions for the differential evolution algorithm
    def predict_fn(xs):
        return predict_classes(xs, x_test[img_id], target_class, 
                               model, target is None)
    
    def callback_fn(x, convergence):
        print(x, target_class, targeted_attack, verbose)
        return attack_success(x, x_test[img_id], target_class, 
                              model, targeted_attack, verbose)
    
    # Call Scipy's Implementation of Differential Evolution
    attack_result = differential_evolution(
        predict_fn, bounds, maxiter=maxiter, popsize=popmul,
        recombination=1, atol=-1, callback=callback_fn, polish=False)

    # Calculate some useful statistics to return from this function
    attack_image = perturb_image(attack_result.x, x_test[img_id])[0]
    prior_probs = predict_one(x_test[img_id])
    predicted_probs = predict_one(attack_image)
    predicted_class = np.argmax(predicted_probs)
    actual_class = y_test[img_id, 0]
    success = predicted_class != actual_class
    cdiff = prior_probs[actual_class] - predicted_probs[actual_class]

    # Show the best attempt at a solution (successful or not)
    helper.plot_image(attack_image, actual_class, class_names, predicted_class)

    return [model.name, pixel_count, img_id, actual_class, predicted_class, success, cdiff, prior_probs, predicted_probs, attack_result.x]

def perturb_image(xs, img):
    # If this function is passed just one perturbation vector,
    # pack it in a list to keep the computation the same
    if xs.ndim < 2:
        xs = np.array([xs])
    
    # Copy the image n == len(xs) times so that we can 
    # create n new perturbed images
    tile = [len(xs)] + [1]*(xs.ndim+1)
    imgs = np.tile(img, tile)
    
    # Make sure to floor the members of xs as int types
    xs = xs.astype(int)
    
    for x,img in zip(xs, imgs):
        # Split x into an array of 5-tuples (perturbation pixels)
        # i.e., [[x,y,r,g,b], ...]
        pixels = np.split(x, len(x) // 5)
        for pixel in pixels:
            # At each pixel's x,y position, assign its rgb value
            x_pos, y_pos, *rgb = pixel
            img[x_pos, y_pos] = rgb
    
    return imgs

def extract_attack(image_id, x, y, r, g, b):
    pixel = np.array([x,y,r,g,b])
    original_image = x_test[image_id]
    image_perturbed = perturb_image(pixel, x_test[image_id])[0]
    true_class = y_test[image_id, 0]
    prior_confidence = predict_one(x_test[image_id])[true_class]
    predicted_probs = predict_one(image_perturbed)
    predicted_class = np.argmax(predicted_probs)
    attacked_xy = np.floor(pixel[0:2]).astype(int)
    attacked_gene = gene_names[ attacked_xy[0]*32 + attacked_xy[1] ]
    attacked_sample = test_sample_names[image_id]
    attacked_expression = r*65536 + g*256 + b
    orig_r, orig_g, orig_b = x_test[image_id][x][y]
    original_expression = orig_r*65536 + orig_g*256 + orig_b
    return(attacked_sample, attacked_gene,  
           original_expression, attacked_expression, 
           true_class, predicted_class,
           original_image, image_perturbed)


f=open(RESULT + "attack_summary","r")
f2=open(RESULT + "attack_summary_extracted","w")

for line in f:
    items=line.split("\t")
    a,b,c,d,e,f,g,h = extract_attack(int(items[0]),int(items[8]),int(items[9]),
                                     int(items[10]),int(items[11]),int(items[12]))
    print(a,b,c,d,e,f, file=f2, sep='\t')
    print(items[0], int(items[0]), a)
    fig, axes = plt.subplots(1, 2, figsize=(8, 6))
    axes[0].axis('off')   
    axes[0].imshow(g, aspect='equal')
    axes[1].axis('off')
    axes[1].imshow(h, aspect='equal')
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0.1)
    plt.margins(5,5)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.savefig(RESULT + "attack_images/" + a + "_" + b + "_orig_attacked.png", bbox_inches = 'tight', pad_inches = 0, dpi=100)
    plt.close()

f2.close()

