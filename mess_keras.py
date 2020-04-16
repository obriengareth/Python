"""

gareth.obrien


"""


##############################################################################
'''
https://pythonprogramming.net/introduction-deep-learning-python-tensorflow-keras/
'''
###############################################################################

import tensorflow.keras as keras
import tensorflow as tf
import matplotlib.pyplot as plt
import numpy as np

print(tf.__version__)
mnist = tf.keras.datasets.mnist
(x_train, y_train),(x_test, y_test) = mnist.load_data()

plt.imshow(x_train[0],cmap=plt.cm.binary)
plt.show()

#It's generally a good idea to "normalize" your data. This typically involves 
#scaling the data to be between 0 and 1, or maybe -1 and positive 1. In 
#our case, each "pixel" is a feature, and each feature currently ranges 
#from 0 to 255. Not quite 0 to 1. Let's change that with a handy utility function:

x_train = tf.keras.utils.normalize(x_train, axis=1)
x_test = tf.keras.utils.normalize(x_test, axis=1)

'''
A sequential model is what you're going to use most of the time. It just means 
things are going to go in direct order. A feed forward model. No going backwards...for now.
Now, we'll pop in layers. Recall our neural network image? Was the input layer flat, 
or was it multi-dimensional? It was flat. So, we need to take this 28x28 image, 
and make it a flat 1x784. There are many ways for us to do this, but keras has a 
Flatten layer built just for us, so we'll use that.
'''


model = tf.keras.models.Sequential()
model.add(tf.keras.layers.Flatten())

# add in 128 hidden layer 1
model.add(tf.keras.layers.Dense(128, activation=tf.nn.relu))
# add in another hidden layer
model.add(tf.keras.layers.Dense(128, activation=tf.nn.relu))
# output layer 10 digits
model.add(tf.keras.layers.Dense(10, activation=tf.nn.softmax))

#compile and build
model.compile(optimizer='adam',loss='sparse_categorical_crossentropy',metrics=['accuracy'])

model.fit(x_train, y_train, epochs=3)


val_loss, val_acc = model.evaluate(x_test, y_test)  # evaluate the out of sample data with model
print(val_loss)  # model's loss (error)
print(val_acc)  # model's accuracy

# save and load model code
#model.save('epic_num_reader.model')
#new_model = tf.keras.models.load_model('epic_num_reader.model')

predictions = model.predict(x_test)

#print(predictions)

print(np.argmax(predictions[5]))
plt.imshow(x_test[5],cmap=plt.cm.binary)
plt.show()

###############################################################################
###############################################################################
###############################################################################

###############################################################################

#https://pythonprogramming.net/convolutional-neural-network-deep-learning-python-tensorflow-keras/?completed=/loading-custom-data-deep-learning-python-tensorflow-keras/

###############################################################################


import tensorflow as tf
from tensorflow.keras.datasets import cifar10
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten
from tensorflow.keras.layers import Conv2D, MaxPooling2D

import pickle

pickle_in = open("X.pickle","rb")
X = pickle.load(pickle_in)

pickle_in = open("y.pickle","rb")
y = pickle.load(pickle_in)

X = X/255.0

model = Sequential()

model.add(Conv2D(256, (3, 3), input_shape=X.shape[1:]))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))

model.add(Conv2D(256, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))

model.add(Flatten())  # this converts our 3D feature maps to 1D feature vectors

model.add(Dense(64))

model.add(Dense(1))
model.add(Activation('sigmoid'))

model.compile(loss='binary_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])

model.fit(X, y, batch_size=32, epochs=3, validation_split=0.3)











