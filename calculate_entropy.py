import math

def get_entropy(profile):   
    total_entropy = [] #will be sum of column totals
    final_entropy = 0
    num = 0
    column = 0
    k = len(profile["A"])
    for j in range (k): #for each column
        for symbol in "ACGT":
            entropy = [] #reset entropy each time
            entropy.append(profile[symbol][j])
            #print(entropy)
            for i in range(len(entropy)):
                #print("entropy " + str(entropy[i]))
                if entropy[i] != 0:
                    total_entropy.append(entropy[i] * math.log(entropy[i], 2))
    for i in range(len(total_entropy)):
        final_entropy += total_entropy[i]      
    return -final_entropy

current_profile = {'A':[0.2,0.2,0,0,0,0,0.9,0.1,0.1,0.1,0.3,0], 'C':[0.1,0.6,0,0,0,0,0,0.4,0.1,0.2,0.4,0.6], 'G':[0,0,1,1,0.9,0.9,0.1,0,0,0,0,0], 'T':[0.7,0.2,0,0,0.1,0.1,0,0.5,0.8,0.7,0.3,0.4]}
print(get_entropy(current_profile))
