from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt
Zeros=0
Zeroz=0
list_of_healthy_tissues=[]
list_of_cancerous_tissues=[]
list_of_healthy_tissues_without_Zeros=[]
list_of_cancerous_tissues_without_Zeros=[]
NamesForHealthy=[]
NamesForCanserous=[]
Correlations_list=[]


# Open file in reading mood

healthy_tissues = open("lusc-rsem-fpkm-tcga_paired.txt","r") 
cancerous_tissues = open("lusc-rsem-fpkm-tcga-t_paired.txt","r") 

# Using for loop to make list of lists by coverting every row to list using .split() then using .append()
# to put them in a list 

for h in healthy_tissues:
    gene = h.split() 
    list_of_healthy_tissues.append(gene)
for c in cancerous_tissues:
    gene = c.split() 
    list_of_cancerous_tissues.append(gene)

# Applying for loop to iterate over the two lists in parallel using zip() and start with the seconed row for 
# each of them to neglect the first row of strings 

for a , b in zip(list_of_healthy_tissues[1:] , list_of_cancerous_tissues[1:]): 
    # Resetting the two counters
    Zeros=0  
    Zeroz=0 
    # Applying for loop ti iterate over the internal lists for every list to make a filteration for rows 
    # which have zeros more than 50% 
    for x , y in zip(a , b):
        if x == '0.0':
            #increament Zeros counter
            Zeros+=1
        if y =='0.0':
            #increament Zeros counter
            Zeroz +=1   
    # Check if numbers of zeros in every list if smaller than 25 to put it in a new list         
    if Zeros <= 25 and Zeroz <= 25 :
        # Adding these lists in new lists
        NamesForHealthy.append(a)
        NamesForCanserous.append(b)
        # Adding these lists in new lists with neglecting the two first columns (for name and id)
        list_of_healthy_tissues_without_Zeros.append(a[2:])
        list_of_cancerous_tissues_without_Zeros.append(b[2:])       

# Converting sublists of strings to sublists of floats

list_of_cancerous_tissues_without_Zeros= [list(map(float, sublist)) for sublist in list_of_cancerous_tissues_without_Zeros] 
list_of_healthy_tissues_without_Zeros= [list(map(float, sublist)) for sublist in list_of_healthy_tissues_without_Zeros]

# Applying for loop to iterate over the new lists to calculate correlation for every gene and 
# putting all correlations in a list  

for a , b in zip(list_of_healthy_tissues_without_Zeros,list_of_cancerous_tissues_without_Zeros):
    correlation,_ = pearsonr(a,b) 
    Correlations_list.append(correlation)

# Getting the names of the genes which have max and min correlation using the index of max and min value 
MaxGeneName=NamesForHealthy[Correlations_list.index(max(Correlations_list))][0]
MinGeneName=NamesForHealthy[Correlations_list.index(min(Correlations_list))][0]  
# Getting the list of healthy tissues which has max correlation healthy tissue 
Xmax=list_of_healthy_tissues_without_Zeros[Correlations_list.index(max(Correlations_list))]
# Getting the list of cancerous tissues which has max correlation healthy tissue 
Ymax=list_of_cancerous_tissues_without_Zeros[Correlations_list.index(max(Correlations_list))]  
# Getting the list of healthy tissues which has min correlation healthy tissue   
Xmin=list_of_healthy_tissues_without_Zeros[Correlations_list.index(min(Correlations_list))]
# Getting the list of cancerous tissues which has min correlation healthy tissue 
Ymin=list_of_cancerous_tissues_without_Zeros[Correlations_list.index(min(Correlations_list))]    
# Plot Positive Correlation
plt.scatter(Xmax, Ymax) 
plt.title(f'A plot to show the Positive correlation between Healthy and Cancerous Tissues \n {MaxGeneName}')
plt.xlabel('Healthy')
plt.ylabel('Cancerous')
plt.plot(np.unique(Xmax), np.poly1d(np.polyfit(Xmax, Ymax, 1))(np.unique(Xmax)), color='yellow')
plt.show()       
# Plot Negative Correlation
plt.scatter(Xmin, Ymin) 
plt.title(f'A plot to show the Nagative correlation between Healthy and Cancerous Tissues \n {MinGeneName}')
plt.xlabel('Healthy')
plt.ylabel('Cancerous')
plt.plot(np.unique(Xmin), np.poly1d(np.polyfit(Xmin, Ymin, 1))(np.unique(Xmin)), color='yellow')
plt.show()       
# Close files
healthy_tissues.close()
cancerous_tissues.close()


