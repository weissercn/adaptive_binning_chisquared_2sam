from __future__ import print_function
import numpy as np
#import matplotlib
#from matplotlib.pylab import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
import time


label_size = 28



################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
mpl.rc('font', family='serif', size=34, serif="Times New Roman")

#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

mpl.rcParams['legend.fontsize'] = "medium"

mpl.rc('savefig', format ="pdf")

mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
mpl.rcParams['figure.figsize']  = 8, 6
mpl.rcParams['lines.linewidth'] = 2

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################




def histo_plot_pvalue(U_0,abins,axlabel,aylabel,atitle,aname):
        bins_probability=np.histogram(U_0,bins=abins)[1]

        #Finding the p values corresponding to 1,2 and 3 sigma significance.
        no_one_std_dev=sum(i < (1-0.6827) for i in U_0)
        no_two_std_dev=sum(i < (1-0.9545) for i in U_0)
        no_three_std_dev=sum(i < (1-0.9973) for i in U_0)

        print(no_one_std_dev,no_two_std_dev,no_three_std_dev)

        with open(aname+"_p_values_1_2_3_std_dev.txt",'w') as p_value_1_2_3_std_dev_file:
                p_value_1_2_3_std_dev_file.write(str(no_one_std_dev)+'\t'+str(no_two_std_dev)+'\t'+str(no_three_std_dev)+'\n')

        #plt.rc('text', usetex=True)
        textstr = '$1\sigma=%i$\n$2\sigma=%i$\n$3\sigma=%i$'%(no_one_std_dev, no_two_std_dev, no_three_std_dev)


        # Making a histogram of the probability predictions of the algorithm. 
        fig_pred_0= plt.figure()
        ax1_pred_0= fig_pred_0.add_axes([0.2,0.15,0.75,0.8]) 
        n0, bins0, patches0 = ax1_pred_0.hist(U_0, bins=bins_probability, facecolor='red', alpha=0.5)
        ax1_pred_0.set_xlabel(axlabel)
        ax1_pred_0.set_ylabel(aylabel)
        #ax1_pred_0.set_title(atitle+ " p values")
        plt.xlim([0,1])

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in upper left in axes coords
        ax1_pred_0.text(0.85, 0.95, textstr, transform=ax1_pred_0.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

        fig_pred_0.savefig(aname+"_p_values_plot.pdf")
        plt.close(fig_pred_0)

def adaptive_binning_2Dplot(bin_boundaries_dict,signed_Scp2_dict,number_of_splits,X_values,title,name,X_min=[0.,0.], X_total_width=[1.,1.]):
	if X_values.shape[1]==2:
		fig = plt.figure()
		ax = fig.add_axes([0.2,0.15,0.75,0.8])
		signed_Scp2_values = []
		for bin_key, bin_Scp2 in signed_Scp2_dict.items(): 
			if str(number_of_splits) in bin_key: 
				signed_Scp2_values.append(bin_Scp2)	
		norm = matplotlib.colors.Normalize(min(signed_Scp2_values),max(signed_Scp2_values))
		my_cmap = cm.get_cmap('Greys')
		
		#ax1.scatter(X_values[:,0], X_values[:,1], c='black',s=5, lw = 0)

		xedges = np.linspace(X_min[0],X_min[0]+X_total_width[0],30)
                yedges = np.linspace(X_min[1],X_min[1]+X_total_width[1],30)
		H, xedges, yedges = np.histogram2d(X_min[0]+X_total_width[0]*X_values[:,0], X_min[1]+X_total_width[1]*X_values[:,1], bins=(xedges, yedges))
		ax.imshow(H, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap='Greys',aspect='auto')


		np.savetxt(name + "_bin_definitions_2D_X_values.txt",X_values)
		boundaries_file = open(name + "_bin_definitions_2D_bin_boundaries_signedScp2.txt", 'w')

		for bin_key, bin_boundary in bin_boundaries_dict.items():
			if str(number_of_splits) in bin_key:
				#print("Boundaries : ",bin_boundary[0,0], bin_boundary[1,0], bin_boundary[0,1],bin_boundary[1,1])
				#cd=ax1.add_patch(Rectangle((bin_boundary[0,0], bin_boundary[1,0]),bin_boundary[0,1]-bin_boundary[0,0] , bin_boundary[1,1]-bin_boundary[1,0],alpha=0.5, edgecolor="magenta",fill=False,linewidth=2))
				#boundaries_file.write("%s\t%s\t%s\t%s\t%s\n" % (bin_boundary[0,0], bin_boundary[1,0], bin_boundary[0,1], bin_boundary[1,1], signed_Scp2_dict[bin_key]))			
				cd=ax.add_patch(Rectangle((X_min[0]+X_total_width[0]*bin_boundary[0,0], X_min[1]+X_total_width[1]*bin_boundary[1,0]),X_min[0]+X_total_width[0]*bin_boundary[0,1]-X_min[0]+X_total_width[0]*bin_boundary[0,0] , X_min[1]+X_total_width[1]*bin_boundary[1,1]-X_min[1]+X_total_width[1]*bin_boundary[1,0],alpha=0.5, edgecolor="magenta",fill=False,linewidth=2))
                                boundaries_file.write("%s\t%s\t%s\t%s\t%s\n" % (X_min[0]+X_total_width[0]*bin_boundary[0,0], X_min[1]+X_total_width[1]*bin_boundary[1,0], X_min[0]+X_total_width[0]*bin_boundary[0,1], X_min[1]+X_total_width[1]*bin_boundary[1,1], signed_Scp2_dict[bin_key]))    

	
		boundaries_file.close()		

		ax.set_xlim(X_min[0],X_min[0]+ X_total_width[0])
		ax.set_ylim(X_min[1],X_min[1]+ X_total_width[1])	
                ax.set_xlabel(r"$x_1$")
                ax.set_ylabel(r"$x_2$")
                #plt.title(name + " bin definitions" )
		fig.savefig(name+"_bin_definitions_2D.pdf")
		plt.close(fig)
		print("The plot {}_bin_definitions_2D.pdf has been made".format(name))
	else:
		print("You can only make this plot in 2D")


def adaptive_binning_2D1Dplot(bin_boundaries_dict,bins_sample01_dict,number_of_splits,X_values,title,name, dim):
        if X_values.shape[1]>1:

		fig = plt.figure()
		ax = fig.add_axes([0.2,0.15,0.75,0.8])

                #ax1.scatter(X_values[:,0], X_values[:,1], c='black',s=5, lw = 0)

                bin_entries0 = []
                bin_entries1 = []
                for bin_key, bin_boundary in bin_boundaries_dict.items():
                        if str(number_of_splits) in bin_key:
                                #print("Boundaries : ",bin_boundary[0,0], bin_boundary[1,0], bin_boundary[0,1],bin_boundary[1,1])
				
				bin_entries0.append(bins_sample01_dict[bin_key][0])
				bin_entries1.append(bins_sample01_dict[bin_key][1])

		bin_middle = np.add(1,range(len(bin_entries0)))
		xwidths = [0.5]*(len(bin_entries0))

		#print("bins_sample01_dict : ", bins_sample01_dict)
                
		#ax.errorbar(bin_middle, bin_entries0, xerr=xwidths, yerr=np.sqrt(bin_entries0), linestyle='', marker='s', markeredgewidth=0.0, markersize=10, color='blue', label=r'$f_1$')
		#ax.errorbar(bin_middle, bin_entries1, xerr=xwidths, yerr=np.sqrt(bin_entries1), linestyle='', marker='o', markeredgewidth=0.0, markersize=10, color='red', label=r'$f_2$')

		ax.errorbar(bin_middle, bin_entries0, yerr=np.sqrt(bin_entries0), linestyle='', marker='s', markeredgewidth=0.0, markersize=10, color='blue', label=r'$f_1$')
                ax.errorbar(bin_middle, bin_entries1, yerr=np.sqrt(bin_entries1), linestyle='', marker='o', markeredgewidth=0.0, markersize=10, color='red', label=r'$f_2$')
		ax.plot((0.,9.),(10000./8.,10000./8.),c="grey",linestyle="--")
		#ax.text(0.05, 0.95,'d = {}'.format(dim), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
		ax.text(0.40, 0.95,'d = {}'.format(dim), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

                ax.set_xlabel("Bin Number")
                ax.set_ylabel("Events")
		#ax = plt.figure().gca()
                #ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		#plt.xticks(range(1,1+len(bin_entries0)))

		ax.set_xlim(0.5,len(bin_entries0)+0.5)
		ax.set_ylim((X_values.shape[0]/(2**number_of_splits))*0.75/2.,(X_values.shape[0]/(2**number_of_splits))*1.25/2.)
		#plt.gca().set_ylim(bottom=0)
		ax.legend(loc='lower right', frameon=False, numpoints=1)
                #plt.title(name + " bin definitions" )
                fig.savefig(name+"_bin_definitions_2D1D.pdf")
                plt.close(fig)
                print("The plot {}_bin_definitions_2D1D.pdf has been made".format(name))
        else:
                print("You can only make this plot in 2 or higher D")



def adaptive_binning_1Dplot(bin_boundaries_dict,data,number_of_splits,title,name):
	if data.shape[1]-1==1:
		data0    = data[data[:,-1]==0]
                data1    = data[data[:,-1]==1]
		#print("data0 : ",data0)
		features_0 = data0[:,:-1]
		features_1 = data1[:,:-1]
		#print("features_0 : ",features_0)

		fig = plt.figure()
		ax = fig.add_axes([0.2,0.15,0.75,0.8])
		#bin_sample01_list = []
		bin_boundaries = []
		bin_boundaries_lower = []
		for bin_key, bin_boundaries_in_dict in bin_boundaries_dict.items():
			if str(number_of_splits) in bin_key:
				#bin_sample01_list.append(bins_sample01)
				bin_boundaries.append(bin_boundaries_in_dict[0][1])
				bin_boundaries_lower.append(bin_boundaries_in_dict[0][0])

		bin_boundaries = sorted(bin_boundaries)
                #bin_boundaries, features = zip(*sorted(zip(bin_boundaries,bin_sample01_list)))
		#print("bin_boundaries : ",bin_boundaries)
		bin_boundaries = [min(bin_boundaries_lower)] + list(bin_boundaries)
		#print("bin_boundaries : ",bin_boundaries)

		#print("bin_sample01_list : ",np.array(bin_sample01_list))

		#bin_sample0 = np.array(bin_sample01_list)[:,0]
                #bin_sample1 = np.array(bin_sample01_list)[:,1]
		#print("features_0.shape : ",features_0.shape,"\nfeatures_1.shape : ",features_1.shape)

		features_0_string, features_1_string, bin_boundaries_string = str(features_0[0][0]), str(features_1[0][0]), str(bin_boundaries[0])

		for i in range(1,features_0.shape[0]):
			features_0_string += "\t"+str(features_0[i][0])			
                for i in range(1,features_1.shape[0]):
                        features_1_string += "\t"+str(features_1[i][0]) 
                for i in range(1,len(bin_boundaries)):
                        bin_boundaries_string += "\t"+str(bin_boundaries[i]) 

		features_boundaries_file = open(name + "_bin_definitions_1D_features_boundaries.txt", 'w')
		features_boundaries_file.write(features_0_string+"\n")
		features_boundaries_file.write(features_1_string+"\n")
		features_boundaries_file.write(bin_boundaries_string+"\n")		
		features_boundaries_file.close()

		#plt.hist(features_0,bins=bin_boundaries,color='red',label='original',alpha=0.4,normed=True)
		#plt.hist(features_1,bins=bin_boundaries,color='blue',label= 'modified',alpha=0.4,normed=True)
                
                hist0, hist0_edges = np.histogram(features_0, bins=bin_boundaries)
                hist1, hist1_edges = np.histogram(features_1, bins=bin_boundaries)
                bin_middle = (hist0_edges[1:] + hist0_edges[:-1]) / 2
		xwidths = (hist0_edges[1:]-hist0_edges[:-1])/2.

		ax.errorbar(bin_middle, hist0, xerr=xwidths, yerr=np.sqrt(hist0), linestyle='', marker='*', markersize=15, color='blue', label=r'$f_1$')
		ax.errorbar(bin_middle, hist1, xerr=xwidths, yerr=np.sqrt(hist1), linestyle='', marker='*', markersize=15, color='red', label=r'$f_2$')

                ax.set_xlim(0.,1.)
                plt.gca().set_ylim(bottom=0)

		ax.set_xlabel(r"$\mathit{x}_{1norm}$")
		ax.set_ylabel("entries")
		ax.legend(loc='best', frameon=False)

                #plt.title(title + " bin definitions" )
		fig.savefig(name+"_bin_definitions_1D.pdf")
		plt.close(fig)
                print("The plot {}_bin_definitions_1D.pdf has been made".format(name))
        else:
                print("You can only make this plot in 1D")
		

def regular_binning_2Dplot(keys, signed_Scp2,single_no_bins,X_values,title,name):
        no_dim = X_values.shape[1]
	if no_dim==2:
		no_bins = [single_no_bins]*no_dim

		fig = plt.figure()
		ax = fig.add_axes([0.2,0.15,0.75,0.8])

                norm = matplotlib.colors.Normalize(np.min(signed_Scp2),np.max(signed_Scp2))
                my_cmap = cm.get_cmap('Greys')

                #ax1.scatter(X_values[:,0], X_values[:,1], c='black',s=5, lw = 0)

                xedges = np.linspace(0.,1.,30)
                yedges = np.linspace(0.,1.,30)
                H, xedges, yedges = np.histogram2d(X_values[:,0], X_values[:,1], bins=(xedges, yedges))
                ax.imshow(H, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap='Greys', aspect='auto')

		np.savetxt(name + "_bin_definitions_2D_X_values.txt",X_values)
                boundaries_file = open(name + "_bin_definitions_2D_bin_boundaries_signedScp2.txt", 'w')

		for key_index, key in enumerate(keys):
			key_list = np.array(key.split(',')).astype(int).tolist()
 

			bin_boundary = np.array([[float(key_list[0])/float(single_no_bins)  ,float(key_list[0]+1)/float(single_no_bins) ],[float(key_list[1])/float(single_no_bins)  ,float(key_list[1]+1)/float(single_no_bins)  ]])
			#print("bin_boundary : ",bin_boundary)
			#print("bin_boundary[0,0] : ",bin_boundary[0,0])
			#print("signed_Scp2 : ", signed_Scp2)


			cd=ax.add_patch(Rectangle((bin_boundary[0,0], bin_boundary[1,0]),bin_boundary[0,1]-bin_boundary[0,0] , bin_boundary[1,1]-bin_boundary[1,0],alpha=0.5, edgecolor="magenta",fill=False,linewidth=2))
	
			boundaries_file.write("%s\t%s\t%s\t%s\t%s\n" % (bin_boundary[0,0], bin_boundary[1,0], bin_boundary[0,1], bin_boundary[1,1], signed_Scp2[key_index]))

		boundaries_file.close()

                ax.set_xlim(0,1)
                ax.set_ylim(0,1)
                ax.set_xlabel(r"$\mathit{x}_{1norm}$")
		ax.set_ylabel(r"$\mathit{x}_{2norm}$")
                #plt.title(title + " bin definitions" )

		fig.savefig(name+"_bin_definitions_2D.pdf")
		plt.close(fig)
                print("The plot {}_bin_definitions_2D.pdf has been made".format(name))
        else:
                print("You can only make this plot in 2D")

def regular_binning_2D1Dplot(keys, bin_points_dict, single_no_bins,X_values,title,name):
        no_dim = X_values.shape[1]
        if no_dim>1:                
                no_bins = [single_no_bins]*no_dim

		fig = plt.figure()
		ax = fig.add_axes([0.2,0.15,0.75,0.8])
                #ax1.scatter(X_values[:,0], X_values[:,1], c='black',s=5, lw = 0)

		bin_entries0 = []
		bin_entries1 = []
                for key_bin_points, bin_points in bin_points_dict.items():
			
			bin_entries0.append(bin_points[0])
			bin_entries1.append(bin_points[1])

		bin_middle = np.add(1,range(len(bin_entries0)))
		xwidths = [0.5]*(len(bin_entries0))		

		print("bin_entries0 : ",bin_entries0)
		print("bin_middle : ",bin_middle)

		ax.errorbar(bin_middle, bin_entries0, xerr=xwidths, yerr=np.sqrt(bin_entries0), linestyle='', marker='*', markersize=15, color='blue', label=r'$f_1$')
		ax.errorbar(bin_middle, bin_entries1, xerr=xwidths, yerr=np.sqrt(bin_entries1), linestyle='', marker='*', markersize=15, color='red', label=r'$f_2$')

                ax.set_xlabel("Bin Number")
                ax.set_ylabel("Entries")


                ax.set_xlim(0.5,len(bin_entries0)+0.5)
                plt.gca().set_ylim(bottom=0)
		#ax = plt.figure().gca()
		#ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		#plt.xticks(range(1,1+len(bin_entries0)))
                #plt.title(title + " bin definitions" )
		ax.legend(loc='best', frameon=False)

                fig.savefig(name+"_bin_definitions_2D1D.pdf")
                plt.close(fig)
                print("The plot {}_bin_definitions_2D1D.pdf has been made".format(name))
        else:
                print("You can only make this plot in 2 or higher D")



def regular_binning_1Dplot(data,single_no_bins,title,name):
        if data.shape[1]-1==1:
                data0    = data[data[:,-1]==0]
                data1    = data[data[:,-1]==1]
                #assert data.shape[0]==10000
		#print("data0 : ",data0)
                features_0 = data0[:,:-1]
                features_1 = data1[:,:-1]

                fig = plt.figure()
		ax = fig.add_axes([0.2,0.15,0.75,0.8])

		bin_boundaries = np.linspace(0,1,single_no_bins+1)

                features_0_string, features_1_string, bin_boundaries_string = str(features_0[0]), str(features_1[0]), str(bin_boundaries[0])

                for i in range(1,features_0.shape[0]):
                        features_0_string += "\t"+str(features_0[i])
                for i in range(1,features_1.shape[0]):
                        features_1_string += "\t"+str(features_1[i])
                for i in range(1,bin_boundaries.shape[0]):
                        bin_boundaries_string += "\t"+str(bin_boundaries[i])

		#assert features_0.shape[0]==10000

                features_boundaries_file = open(name + "_bin_definitions_1D_features_boundaries.txt", 'w')
                features_boundaries_file.write(features_0_string+"\n")          
                features_boundaries_file.write(features_1_string+"\n")
                features_boundaries_file.write(bin_boundaries_string+"\n")
		features_boundaries_file.write(str(features_0.shape[0])+"\t"+str(features_1.shape[0]))
		features_boundaries_file.close()
	
                #plt.hist(features_0,bins=bin_boundaries,color='red',label='original',alpha=0.4,normed=True)
                #plt.hist(features_1,bins=bin_boundaries,color='blue',label= 'modified',alpha=0.4,normed=True)
	
		hist0, hist0_edges = np.histogram(features_0, bins=bin_boundaries)
		hist1, hist1_edges = np.histogram(features_1, bins=bin_boundaries)
		bin_middle = (hist0_edges[1:] + hist0_edges[:-1]) / 2
		xwidths = (hist0_edges[1:]-hist0_edges[:-1])/2.

		hist0_sqrt = np.sqrt(hist0)
		hist1_sqrt = np.sqrt(hist1)

		scaling0 = float(xwidths.shape[0])/sum(hist0)
		scaling1 = float(xwidths.shape[0])/sum(hist1)

		hist0 = np.multiply(hist0,scaling0)
		hist1 = np.multiply(hist1,scaling1)
                hist0_sqrt = np.multiply(hist0_sqrt,scaling0)
                hist1_sqrt = np.multiply(hist1_sqrt,scaling1)

		ax.plot((0.,1.),(1.,1.),c="grey",linestyle="--")	

		ax.errorbar(bin_middle, hist0, xerr=xwidths, yerr=hist0_sqrt, linestyle='', marker='*', markersize=15, color='blue', label=r'$f_1$')
		ax.errorbar(bin_middle, hist1, xerr=xwidths, yerr=hist1_sqrt, linestyle='', marker='*', markersize=15, color='red', label=r'$f_2$')

                ax.set_xlim(0.,1.)
                plt.gca().set_ylim(bottom=0)
		ax.set_xlabel(r"$\mathit{x}_{1norm}$")
		ax.set_ylabel("entries")
		#plt.title(title + " bin definitions" )
		ax.legend(loc='lower right', frameon=False, numpoints=1)

                fig.savefig(name+"_bin_definitions_1D.pdf")
		plt.close(fig)
                print("The plot {}_bin_definitions_1D.pdf has been made".format(name))
        else:
                print("You can only make this plot in 1D")


























