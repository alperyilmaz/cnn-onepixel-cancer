import pickle
import numpy as np
import pandas as pd
import matplotlib
from keras.datasets import cifar10
from keras import backend as K
import tensorflow as tf
import keras
from datetime import datetime
from differential_evolution import differential_evolution
import helper
import random

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


### attack functions ###

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
    # max gene expression was found to be 5779850, corresponding to 88 in red
    # however, when we count how many times each R value is seen, there's no need to go over 2
    # gunzip -c TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz  | _filtergenes |  _pivotlonger | _filtersamples | awk '{printf"%s\t%s\t%.0f\n",$1,$2,$3}' | awk '{r=and(rshift($3,16),255); count[r]++;}END{for(val in count){print val,count[val]}}'
    # 0 18266327
    # 1 35887
    # 2 6472
    # 3 2292
    # 4 932
    # 5 507
    # 6 262
    # 7 177
    # 8 98
    bounds = [(0,32), (0,32), (0,2), (0,256), (0,256)] * pixel_count
    
    # Population multiplier, in terms of the size of the perturbation vector x
    popmul = max(1, popsize // len(bounds))
    
    # Format the predict/callback functions for the differential evolution algorithm
    def predict_fn(xs):
        return predict_classes(xs, x_test[img_id], target_class, 
                               model, target is None)
    
    def callback_fn(x, convergence):
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
    #helper.plot_image(attack_image, actual_class, class_names, predicted_class)

    return [model.name, pixel_count, img_id, actual_class, predicted_class, success, cdiff, prior_probs, predicted_probs, attack_result.x]
    


def attack_all(models, samples=500, pixels=(1), targeted=False, 
               maxiter=75, popsize=400, verbose=False):
    results = []
    for model in models:
        model_results = []
        valid_imgs = correct_imgs[correct_imgs.name == model.name].img
        img_samples = np.random.choice(valid_imgs, samples, replace=False)
        
        for pixel_count in pixels:
            for i, img_id in enumerate(img_samples):
                print('\n', model.name, '- image', img_id, '-', i+1, '/', len(img_samples))
                targets = [None] if not targeted else range(2)
                
                for target in targets:
                    if targeted:
                        print('Attacking with target', class_names[target])
                        if target == y_test[img_id, 0]:
                            continue
                    result = attack(img_id, model, target, pixel_count, 
                                    maxiter=maxiter, popsize=popsize, 
                                    verbose=verbose)
                    model_results.append(result)
                    
        results += model_results
        helper.checkpoint(results, targeted)
    return results
    
### ###

RESULT="results/"
MODEL="model/"

seed_no=random.randint(1,10000)
np.random.seed(seed_no)

class_names = ["Normal","Tumor"]
data = np.load(RESULT + "tcga_gtex_genes_data.npz")
x_train = data["x_train"]
x_test = data["x_test"]
y_train = data["y_train"].reshape(-1,1)
y_test = data["y_test"].reshape(-1,1)

print("[Reading model..]")
tcganet = keras.models.load_model(MODEL + "tcgamodel.h5")
models = [tcganet]

network_stats, correct_imgs = helper.evaluate_models(models, x_test, y_test)
correct_imgs = pd.DataFrame(correct_imgs, columns=['name', 'img', 'label', 'confidence', 'pred'])
network_stats = pd.DataFrame(network_stats, columns=['name', 'accuracy', 'param_count'])

network_stats.to_csv(RESULT + "model_stats.csv")

print("[Starting the attack..]")
now = datetime.now()
dt_string = now.strftime("%Y%m%d-%H%M%S")

untargeted = attack_all(models, samples=800, targeted=False, pixels=[1])

columns = ['model', 'pixels', 'image', 'true', 'predicted', 'success', 'cdiff', 'prior_probs', 'predicted_probs', 'perturbation']

results_table = pd.DataFrame(untargeted, columns=columns)

results_table.to_csv(RESULT + "attack_results/" + "attack_results_seed" + str(seed_no) + "_" + dt_string + ".csv")
