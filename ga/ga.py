# -*- coding: utf-8 -*-
"""
Genetic algarithum for minimum problem!
"""
import itertools
import math
import random

import Levenshtein as le
import pandas as pd
import progressbar
from matplotlib import pyplot as plt


# target function input
class fitness_function:
    def __init__(self):
        self.lower_edge = 1
        self.upper_edge = 1000

    def func(self, x): # for minimum +, for maximum -
        return math.sin(x)+5*math.cos(x)-math.exp(-0.1*x)

# ga main function
class ga(fitness_function):
    def __init__(self):
        super().__init__()
        self.dna_length = 20  # time of 2
        self.population = 100  # numbers of DNA in one generation chrom
        self.generation = 200  # times to loop for best chrom
        if self.population % 2 != 0:
            self.population += 1  # population shall be dual
        self.elite_fitness_num = round(self.population/10) # numbers of elite to store base on fitness
        self.elite_expectation_num = round(self.population/10) # numbers of elite to store base on expectation
        self.relative_valve = 4/self.dna_length  # factor to multiply relative counts, fitness - relative_valve * relative counts
        self.mutation_possibility = 1/self.population  # chance of mutation for each dna unit
        self.exit_flag = self.generation//3  # if generation//3 numbers of result are same, then exit the loop
        self.main_process()  # main running process
        #___________________________________________________________________

    def main_process(self):
        bar = progressbar.ProgressBar(0, self.generation)
        self.first_chrom()  # generate the first chrom
        self.chrom = self.fchrom  # shift the first chrom to common various chrom
        self.result = []  # recording result of function
        for i in range(0, self.generation):
            self.fitness()  # calculate the fitness
            self.result.append(self.chrom.ix[0, 'RESULT'])  # store the best fitness result
            # fitness storage into the chrom dataframe in cloumn FITNESS
            self.elite_storage_fitness()  # storage numbers of elite with best fitness, numbers according to self.elite_fitness_num
            # store into self.elite_fitness
            self.relative()  # count the DNV similar condition with other
            self.expectation()  # calculate the expectation
            self.elite_storage_expectation()  # storage the elite expectation
            # store into self.elite_expectation
            self.chrom = self.new_generation()  # create new generation
            bar.update(i)
            # exit loop if all result same within numbers of self.exit_flag
            if i > self.exit_flag:
                if self.result[i]==self.result[i-self.exit_flag]:
                    print()  #print enter
                    print('Result fix in 1/3 of generation, exit loop!')
                    break   # break loop of generation         
        self.fitness()  # calculate the final generation fitness, with RESULT and FITNESS, sorted by fitness
        self.result.append(self.chrom.ix[0, 'RESULT']) # append the final best result after loop
        best_x = self.mapping(int(self.chrom.ix[0,'DNA'],2))
        plt.plot(list(range(len(self.result))),self.result)  #plot the result in generation
        plt.show()
        print()  #print enter
        print('When x = {x}, Minimum result = {y}'.format(x=best_x, y=self.result[::-1][0]))
     
    
    def first_chrom(self):
        # create first chrom and turn to binary in pandas form: dec to bin, delete '0b', fill zero to DNA length, also create the pandas form
        self.fchrom = pd.DataFrame([str(bin(random.randint(0, 2**self.dna_length)).replace('0b', '')).zfill(self.dna_length) for i in range(self.population)],columns=['DNA'])

    def mapping(self,x):  # to map the DNA range to fitness function range of values
        x_value = (self.upper_edge-self.lower_edge)/2**self.dna_length*x+self.lower_edge
        return x_value

    def mapminmax(self, x):  # mapminnmax, turn list from 0 to 1; x should be a list

        if min(x) == max(x):
            y = list(itertools.repeat(1,len(x)))
        else:
            y = []
            for i in x:
                y.append(1-(i-min(x))/(max(x)-min(x))) # less result better fitness
        return y

    def fitness(self):  # calculate result and fitness of chrom, and concat to the pandas of chrom
        result = []  # recording result of function
        for i in self.chrom.DNA:
            x_value = self.mapping(int(i, 2))  # convert bin to dec and mapping to the value range of target function
            result.append(self.func(x_value))  # append the result to list
        fitness = self.mapminmax(result)  # mapminmax result to fitness, max to 1, min to 0
        fitness = pd.DataFrame([result, fitness], index=['RESULT', 'FITNESS'])  #result and fitness to pandas
        fitness = fitness.T  # row to columns, in order to concat
        self.chrom = pd.concat([self.chrom, fitness], axis=1, join='outer')  # concat the result right behind the chrom
        self.chrom = self.chrom.sort_values(["FITNESS"], ascending=False)  # ascending the chrom by fitness from upper to lower
        self.chrom.index = list(range(len(self.chrom.index)))  # reindex after sorting

    def elite_storage_fitness(self):
        # here terms in fitness function already sorted the fitness
        self.elite_fitness = self.chrom.DNA.head(self.elite_fitness_num)  # store number of top fitness chrom

    def relative(self):
        relative_factor = []  # empty the recording list for all DNV in chrom
        for dna in self.chrom.DNA:  # loop, from 1st DNV to last
            relative_num = 0 # similar count recording for one single DNV
            for cmp_dna in self.chrom.DNA:  # loop DNV for every other DNV, for counting similar relative
                relative_num += le.distance(dna, cmp_dna)  # Levenshtein distance, for counting similar of two strings
            relative_factor.append(relative_num)  # append the counting result for one single DNV
        relative_factor = list(map(lambda x: 1 - x/max(relative_factor), relative_factor))  # similar to 1, differ to 0
        relative_factor = pd.DataFrame(relative_factor, columns=['RELATIVE'])  # change to pandas
        self.chrom = pd.concat([self.chrom, relative_factor], axis=1, join='outer')  # concat to chrom

    def expectation(self):  # calculate the expectation, in order to distinguish a better fitness DNA whether with many relative, if so than lower the possibility
        expectation = self.chrom.FITNESS-self.chrom.RELATIVE*self.relative_valve  # as the script
        expectation = pd.DataFrame(expectation, columns=['EXPECTATION'])  # turn to pandas
        self.chrom = pd.concat([self.chrom, expectation], axis=1, join='outer')  # join to chrom

    def elite_storage_expectation(self): # storage the top numbers of elite of expectation
        self.chrom = self.chrom.sort_values(["EXPECTATION"], ascending=False)  # ascending the chrom by expectation from upper to lower
        self.chrom.index = list(range(len(self.chrom.index)))  # reindex after sorting
        self.elite_expectation = self.chrom.DNA.head(self.elite_expectation_num)  # store number of top fitness chrom

    def wheel_gambling(self): # wheel gambling, to locate one lucky DNA which prepare to marry :)
        possibility = list(itertools.accumulate(self.chrom.EXPECTATION))  # possibility index of accumulation in expectation
        gamble = random.uniform(0, max(possibility))  # random float number between 0 to max accumulation of possibility
        pos_pin = 0 # a pin to locate the random number drop to which single dna
        for i in possibility:  # search the accumulated expectation
            if gamble <= i:  # when if the random point lower than the accumulate expectation, it means it felt to the single DNA
                break # then break the searching
            pos_pin += 1 # if not, then pin to next one
        return self.chrom.ix[pos_pin, 'DNA']  # return the lucky one, with only DNA columns

    def random_dna_section(self):  #pick up random section of dna
        exitflag = True
        while exitflag:  # continue to generate the cross range until they not same
            looks_like = [random.randint(0, self.dna_length+1) for i in range(2)]  # cross range, +1 due to index to include the last one in DNA
            looks_like.sort() #make sure the small one in fore
            if looks_like[0] != looks_like[1]:  # two cross point are not same
                exitflag = False  # then exit loop
        return looks_like

    def coitus(self, x):  # cross the lucky dna, dna string in list to x
        looks_like = self.random_dna_section()  # random section of dna pin
        x_section = []  # for recording dna cross section take out, list for both
        for i in x:
            x_section.append(i[looks_like[0]:looks_like[1]])  # dna cross section take out, list for both
        x[0] = x[0].replace(x[0][looks_like[0]:looks_like[1]], x_section[1])  # cross, the first dna section replace by section of the second dna
        x[1] = x[1].replace(x[1][looks_like[0]:looks_like[1]], x_section[0])  # vice versa
        return x

    def marriage(self):  # pickup two lucky dna and generate two kid crossed
        lucky_dna = [self.wheel_gambling() for i in range(2)]  # pickup two lucky dna
        baby_dna = self.coitus(lucky_dna)  # get two new generation
        return baby_dna

    def new_generation(self):  # create new generation
        loop_generation = int((self.population - self.elite_fitness_num - self.elite_expectation_num)/2)  # time for loop to create new chrom other than the elite
        new_chrom = []  #recorder for new chrom
        for i in range(loop_generation):
            new_generation = self.marriage()  # create two cross new chrom
            for j in new_generation:  # assign two to new chrom one by one
                new_chrom.append(j)
        new_chrom = pd.DataFrame(new_chrom)  # new cross chrom for pandas
        new_chrom = pd.concat([new_chrom, self.elite_fitness, self.elite_expectation])  # concat new chrom, and both elite
        new_chrom.columns = ['DNA']  # rename the DNA
        new_chrom.index = list(range(len(new_chrom)))  # reindex the chrom pandas
        new_chrom = self.mutation(new_chrom)
        return new_chrom

    def mutation(self, x):  # mutation, to keep the chrom variable, input x to be chrom pandas
        for i in range(len(x)):  # loop all DNA
            if random.random() <= self.mutation_possibility: # if random less than the possibility of mutation, than mutation
                looks_like = self.random_dna_section()  # mutation section of dna
                mutate_dna = "".join([random.choice(['0', '1']) for i in range(looks_like[1]-looks_like[0])])  # loop for generate mutate section dna
                x.ix[i, 'DNA'] = x.ix[i, 'DNA'].replace(x.ix[i, 'DNA'][looks_like[0]:looks_like[1]], mutate_dna)
        return x


myfunc = ga()


