import coursework2 as module
import os
import shutil

# Part one of the assignment
print('Generating file with alpha-carbon atom entries.')
module.atoms("5kkk.pdb", atomtype = 'atom') # Creates a file named 'atom_file' with just the CA atom entries

print('Running parsePDB() ...')  # Calls parsePDB and outputs result to console
module.parsePDB()
print('Writing part one results to answer file.')

# Part two of the assignment
print('\nRunning distance_to_centroid() ...')
module.distance_to_centroid("atomfile.txt")
print('Writing part two results to answer file.')

# Part three of the assignment
module.fe_neighbours()
print('\nCalling search_heterogen() to find find the coordinates of the Fe atom... ')
print('Running atom_neighbours() to find the 5 nearest neighbours to the Fe atom...')
print('Writing part three results to answer file.')

# Cleaning up the files in the current directory
print('\nCleaning up intermediary text files into \"TextFiles\" folder.')
temp_destination = os.path.dirname(os.path.abspath(__file__)) + "/TextFiles"
current_dir = os.path.dirname(os.path.abspath(__file__))
files = ['atomfile.txt', 'atomheterofile.txt', 'heterofile.txt', 'phil_atomfile.txt',
         'phob_atomfile.txt', 'phob_atomheterofile.txt', 'phob_heterofile.txt',
         'phil_atomheterofile.txt', 'phil_heterofile.txt']

for file in files:
    source = current_dir + "/" + file
    destination = temp_destination + "/" + file
    shutil.move(source, destination)

# Moving the generated plots to a separate folder
print('Saving generated plots to \"Plots\" folder.')
temp_destination = os.path.dirname(os.path.abspath(__file__)) + "/Plots"
current_dir = os.path.dirname(os.path.abspath(__file__))
files = ['histograms.pdf', 'boxplot.pdf', '3Dscatter.pdf', 'xy2Dscatter.pdf',
         'xz2Dscatter.pdf', 'yz2Dscatter.pdf', 'xy_phobic_KDE.png',
         'xz_phobic_KDE.png', 'yz_phobic_KDE.png', 'xy_philic_KDE.png',
         'xz_philic_KDE.png', 'yz_philic_KDE.png', 'plot_neighbours.pdf']

for file in files:
    source = current_dir + "/" + file
    destination = temp_destination + "/" + file
    shutil.move(source, destination)
