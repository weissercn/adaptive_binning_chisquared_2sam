#!/usr/bin/python
from chi2_regular_binning import *
from chi2_adaptive_binning import *
from chi2_plots import *
import os

def example_call():
	orig_name="gauss__1_0__0_9"
	title="Gauss 0.1 0.75"
	dim_list = [2,3]
	systematics_fraction = 0.01

	comp_file_list_list = []
	adaptive_binning=True
	for dim_data in dim_list:
		comp_file_list=[]
		for i in range(1):
			comp_file_list.append((os.environ['learningml']+"/GoF/data/gaussian_same_projection_on_each_axis/gauss_data/gaussian_same_projection_on_each_axis_redefined_{0}D_10000_0.0_1.0_1.0_{1}.txt".format(dim_data,i),os.environ['learningml']+"/GoF/data/gaussian_same_projection_on_each_axis/gauss_data/gaussian_same_projection_on_each_axis_redefined_{0}D_10000_0.0_1.0_0.9_{1}.txt".format(dim_data,i)))
		comp_file_list_list.append(comp_file_list)
	if adaptive_binning==True:
		number_of_splits_list = [1,2,3]
		#chi2_adaptive_binning.chi2_adaptive_binning(orig_name, dim_list, comp_file_list_list,number_of_splits_list)
		chi2_adaptive_binning_wrapper(title,orig_name, dim_list, comp_file_list_list,number_of_splits_list,systematics_fraction)
	else:
		single_no_bins_list=[2,3,5]
		chi2_regular_binning_wrapper(title,orig_name, dim_list, comp_file_list_list,single_no_bins_list,systematics_fraction)	

if __name__ == "__main__":
	example_call()
