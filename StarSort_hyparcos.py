# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import random

def some(x,n):
    return x.ix[random.sample(x.index, n)]

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def classifySpectrum(spectrum):
    if 'VII' in spectrum or 'D' in spectrum:
        return ['whiteDwarf', 'purple']
    elif 'VI' in spectrum or 'sd' in spectrum:
        return ['subDwarf' , 'grey']
    elif 'IV' in spectrum:
        return ['subGiant' , 'green']
    elif 'III' in spectrum:
        return ['giant' , 'red']
    elif 'II' in spectrum:
        return ['brightGiant' , 'blue']
    elif 'Ib' in spectrum or 'Iab' in spectrum or 'Ia' in spectrum:
        return ['superGiant' , 'black']
    elif 'V' in spectrum:
        return ['mainSequence' , 'orange']
    else:
        return ['NA','NA']
    
blueIndex = -0.33
redIndex  = 1.64
onlyVisibles = True
subSample = False
sampleSize = 6667 #6667

plt.close('all')
df = pd.read_table('hip_main.dat', header = None, skiprows = 0, sep = '|', usecols = [1,11,32,34,40,44,76], names = ['Hindex','paralax', 'BTmag','VTmag', 'ci', 'hipMag', 'spect'])
df  = df[df.spect.str.strip().str.len() != 0]

df['mag'] = df['hipMag']
df = df[df['mag'].apply(is_number)] 
df = df[df['ci'].apply(is_number)] 
df['ci'] = df['ci'].apply(float)
df = df[df['ci'].apply(lambda x: x >= blueIndex)] 
df = df[df['ci'].apply(lambda x: x <= redIndex)] 

df['mag'] = df['mag'].apply(float)

df = df[df['paralax'].apply(is_number)]   
df = df[df['paralax'].apply(lambda x: float(x) > 0)]
df['dist'] = df['paralax'].apply(lambda x: 1000/float(x))
df['absmag'] = df.apply(lambda row:  (row['mag']) - 5*(np.log(row['dist']/10)), axis=1)


print('total stars: ' + str(len(df)))


df['starType'] = df['spect'].apply(lambda x: classifySpectrum(str(x)[2:])[0])
df['plotColor'] = df['spect'].apply(lambda x: classifySpectrum(str(x)[2:])[1])

df = df[df['starType']!= 'NA']
if onlyVisibles:
    df = df[df['mag'] <= 6.5]
    print('total visible stars: ' + str(len(df)))
if subSample:
    df = df.sample(sampleSize)
    
print('number stars in sample: ' + str(len(df)))
#add the sun
df = df.append(pd.DataFrame([[-99]+(len(df.columns)-1)*[0]] ,  columns = list(df.columns)))
df[df['Hindex']==-99]['dist'] = 9999
df = df.set_index('Hindex')
df.set_value(-99,'ci', 0.631)
df.set_value(-99,'absmag', 4.83)
df.set_value(-99,'spect', 'xxV')
df.set_value(-99,'mag', -3)
df.set_value(-99,'starType',classifySpectrum('V')[0])
df.set_value(-99,'plotColor',classifySpectrum('V')[1])

#remove star types
#df = df[df['starType'] != 'whiteDwarf'] 
#df = df[df['starType'] != 'giant'] 
#df = df[df['starType'] != 'brightGiant'] 
#df = df[df['starType'] != 'subDwarf'] 
#df = df[df['starType'] != 'superGiant'] 
#df = df[df['starType'] != 'subGiant'] 
#df = df[df['starType'] != 'main'] 

#df = df[df['ci'] > 1.5] 


colors = df['ci'].tolist()
dists  = df['dist'].tolist()
spectrums = df['spect'].tolist() 
absMags = df['absmag'].tolist() 
starTypes = df['starType'].tolist()
plotColors = df['plotColor'].tolist()



fig = plt.figure(2)
plt.clf()
ax = fig.add_subplot(111)
for ind, i in enumerate(np.unique(starTypes)):
    print(i)
    curPlotColors = [plotColors[ind]  for ind, j in enumerate(starTypes) if i == j]
    curColors = [colors[ind]  for ind, j in enumerate(starTypes) if i == j]
    curAbsMags   = [absMags[ind]  for ind, j in enumerate(starTypes) if i == j]
    plt.scatter(curColors, curAbsMags, s = 7, label = i, color = curPlotColors[0]) 
#add sol
plt.scatter(df.loc[-99]['ci'], df.loc[-99]['absmag'], s = 300, label = 'sol', color = 'red', marker = '.') 
plt.xlabel('color index')
plt.ylabel('abs Mag')
plt.legend()
if onlyVisibles:
    plt.title('H-R diagram of naked eye stars')
elif subSample:
    plt.title('H-R diagram of subsample of hipparcos stars')
else:
    plt.title('H-R diagram of all hipparcos stars')
    
ax.annotate('sol', xy= (df.loc[-99]['ci'], df.loc[-99]['absmag']), xytext = (df.loc[-99]['ci']-0.2, df.loc[-99]['absmag']+5),
            arrowprops=dict( facecolor='black', shrink=0.1)
            )
plt.axvline(x=blueIndex, c = 'blue', linestyle = 'dashed')
plt.axvline(x=redIndex, c = 'red', linestyle = 'dashed')
ax = plt.gca()

extraTicks=[blueIndex, redIndex]
extraTickLabels = ['blue', 'red']
ax.set_xticks( list(ax.get_xticks()) + extraTicks)

fig.canvas.draw()
labels = [item.get_text() for item in ax.get_xticklabels()]

labels[-2] = 'blue'
labels[-1] = 'red'
ax.set_xticklabels( labels )

plt.gca().invert_yaxis()

starCounts = Counter(starTypes)
totalStars = sum(starCounts.values())
starPercentage = [100 * float(i) / totalStars for i in starCounts.values()]

fig = plt.figure(3)
plt.clf()
ax = fig.add_subplot(111)
#fig.canvas.set_window_title('Figure 4') 
ax.bar(range(len(starCounts)), starPercentage, align='center')
#ax.bar(range(len(starCounts)), starPercentage)

ax.set_xticks(range(len(starCounts)))
ax.set_xticklabels(list(starCounts.keys()))
ax.set_ylabel('percent of total')
if onlyVisibles:
    ax.set_title('histogram of visible stars')
elif subSample:
    ax.set_title('histogram of subsample of hipparcos stars')
else:
    ax.set_title('histogram of all Hipparcos stars')
ax.set_ylim(ymax = max(starPercentage) + 5)
rects = ax.patches
for ind, rect in enumerate(rects):
    ax.text(rect.get_x() + rect.get_width()/2, rect.get_height()+ 0.3 ,str(round(starPercentage[ind],2)) + ' %'   )
if onlyVisibles:
    plt.savefig("AllHipparcos_histogram.png", dpi=600)
else:
    plt.savefig("Visibles_histogram.png", dpi=600)
    
    
fig = plt.figure(4)
fig.clf()
ax = fig.add_subplot(111)
maxDist = max(dists)
sigDist = np.std(dists)

plt.set_cmap('rainbow')
plot = ax.scatter(colors, absMags, s = 100, label = i, c = [np.log(i) if i > 0 else 0 for i in dists], marker = '.') 
cax = fig.add_axes([0.91, 0.2, 0.05, 0.5])
#plt.gray()
#ax.set_cmap('Spectral')
fig.colorbar(plot, cax=cax, orientation='vertical')
ax.invert_yaxis()
ax.set_xlabel('color index')
ax.set_ylabel('abs Mag')
ax.annotate('sol', xy= (df.loc[-99]['ci'], df.loc[-99]['absmag']), xytext = (df.loc[-99]['ci']-0.2, df.loc[-99]['absmag']+5),
            arrowprops=dict( facecolor='black', shrink=0.1)
            )
if onlyVisibles:
    ax.set_title('distance of visible stars')
elif subSample:
    ax.set_title('distance of subsample of hipparcos stars')
else:
    ax.set_title('distance of all Hipparcos stars')
    
    
fig = plt.figure(5)
plt.hist([np.log(i) if i > 0 else 0 for i in dists], 50, density=1, facecolor='blue', alpha=0.75)
plt.xlabel('log distance in parsecs')
distHist = np.histogram([np.log(i) if i > 0 else 0 for i in dists], bins = 50)
plt.title('histogram of distance\n most common:' + str(round(3.26*np.exp(distHist[1][np.argmax(distHist[0])]),2)) + ' light years')
