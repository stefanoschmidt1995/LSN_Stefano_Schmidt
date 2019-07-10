#*************************************************************************
#STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
#EXERCISE 11 - HELPER CODE 
#*************************************************************************
import numpy as np
import keras
import tensorflow
import matplotlib.pyplot as plt

from keras.models import Sequential
from keras.layers import Dense, Activation
from keras import backend as K
from keras.utils import get_custom_objects

def create_training_set(N_train, N_valid, sigma, function = None, seed =0, ndim=1, x_range =1):
	"""
	Function to create a dataset out of function (default = straight line).
	Input:	N_train = float
		sigma (float) = noise for observation
		function (lambda type) = scalar function to build the observation from ((n,j) -> (n,))
		ndim (int) = number of dimensions for data (default = 1)
		x_range (float) = float value s.t. x_train is in [-x_range,x_range]
	Output: some np.array
		x_train, y_train, x_val, y_valid, y_true
			y_true = f(x_valid)
	"""
	np.random.seed(seed)
	if function is None:
		m = 2 # slope
		b = 1 # intersect
		function = lambda x: m*x+b

	x_train = np.random.uniform(-x_range, x_range, (N_train,ndim))
	x_valid = np.random.uniform(-x_range, x_range, (N_valid,ndim))
	
	y_train = np.zeros((N_train,))
	y_true = np.zeros((N_valid,))
	y_valid = np.zeros((N_valid,))

	for i in range(N_train):
		y_train[i] = np.random.normal(function(x_train[i,:]), sigma)

	for i in range(N_valid):
		y_true[i] = function(x_valid[i,:])
		y_valid[i] = np.random.normal(function(x_valid[i,:]), sigma)
	return x_train, y_train, x_valid, y_valid, y_true


def fit_NN_model(x_train, y_train, N_epochs, model=None ,valid_data = None):
	"""
	Fit x_train, y_train with a keras model. The model is fitted with sgd. This function is a wrapper to model.compile and model.fit from keras.
	Input:	x_train, y_train (np.array (n,k) (n,))
		N_epochs (float)
		model (keras_model) = keras model to fit
		valid_data (tuple, (x_valid, y_valid)) = couple with a validation set to monitor the performances while training
	Output: list of objects [model, history]
		model = fitted keras model 
		history = keras history object returned by the fitting procedure
	"""
	#creating NN to fit
	input_size = 1
	if x_train.ndim >1:
		input_size= x_train.shape[1]

	if model is None:
		model = keras.Sequential()
		model.add(Dense(1, input_shape=(input_size,)))

	# compile the model choosing optimizer, loss and metrics objects & fitting
#	opt = keras.optimizers.SGD(lr=0.01, momentum=0.01, decay=0.1, nesterov=False)
	opt = 'sgd'
	model.compile(optimizer=opt, loss='mse', metrics=['mse'])
	history = model.fit(x=x_train[:,0:2], y=y_train, batch_size=32, epochs=N_epochs, shuffle=True, verbose=0, validation_data=valid_data)
	return [model,history]


def plot_predictions(model, x_test, y_true, x_train=None, y_train=None):
	"""
	Plot some predictions of a given model and compare it with the true model.
	Input:	model (keras model type)=model that is fitted
		x_test, y_true (np.array) = some points in a test set along with their TRUE y values
		x_train, y_train (np.array) = can be also plotted  the training set used for fitting the model
	Output: None
	"""
	loss = model.evaluate(x_test, y_true, batch_size=32, verbose=0)[0]
	x_predicted = np.random.uniform(-1, 1, 100)
	y_predicted = model.predict(x_predicted)
	plt.scatter(x_predicted, y_predicted, color='r', label='predicted')
	if x_train is not None and y_train is not None:
		plt.scatter(x_train, y_train, color='b', label='training')
	plt.scatter(x_test, y_true, label = 'true')
	plt.title("Predicted vs. true values (loss = "+str(np.round(loss,3))+")")
	plt.grid(True)
	plt.legend()
	plt.show()
	return

def plot_3D_predictions(model, x_test, y_true, x_train=None, y_train=None):
	"""
	Plot some predictions of a given model with 2 dim input and compare it with the true model.
	Input:	model (keras model type)=model that is fitted
		x_test, y_true (np.array) = some points in a test set along with their TRUE y values
		x_train, y_train (np.array) = can be also plotted  the training set used for fitting the model
	Output: None
	"""
	from mpl_toolkits.mplot3d import Axes3D
	fig = plt.figure()
	ax = Axes3D(fig)
	ax.set_xlabel('x')
	ax.set_ylabel('y')

	loss = model.evaluate(x_test, y_true, batch_size=32, verbose=0)[0]
	x_predicted = np.random.uniform(-1.5, 1.5, size=(100,2))
	y_predicted = model.predict(x_predicted)
	ax.scatter(x_predicted[:,0],x_predicted[:,1], y_predicted, color='r', label='predicted')
	if x_train is not None and y_train is not None:
		ax.scatter(x_train[:,0],x_train[:,1], y_train, color='b', label='training')
	ax.scatter(x_test[:,0], x_test[:,1], y_true, color='g', label='true')
	plt.title("Predicted vs. true values (loss = "+str(np.round(loss,3))+")")
	ax.legend();
	ax.view_init(10, 30)
	plt.show()
	return







