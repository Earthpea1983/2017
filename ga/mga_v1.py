# -*- coding: utf-8 -*-
"""
Genetic algarithum for minimum problem, with mutipul various!
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
        self.lower_edge = [-1000,-1000,-1000]
        self.upper_edge = [1000,1000,1000]
        self.various_num = len(self.lower_edge)

    def func(self, x): # for minimum +, for maximum -, x shall be a list
        return x[0]**2+x[1]**2-4*x[2]**2+x[2]*math.sin(x[2])

# ga main function
class ga(fitness_function):
    def __init__(self):
        super().__init__()
        self.dna_length = 22  # time of 2
        self.population = 300  # numbers of DNA in one generation chrom
        self.generation = 300  # times to loop for best chrom
        if self.population % 2 != 0:
            self.population += 1  # population shall be dual
        self.elite_fitness_num = round(self.population/10)  # numbers of elite to store base on fitness
        self.elite_expectation_num = round(self.population/10)  # numbers of elite to store base on expectation
        self.relative_valve = 1/3  # factor to multiply relative counts, fitness - relative_valve * relative counts
        self.mutation_possibility = 1/self.population  # chance of mutation for each dna unit
        self.exit_flag = self.generation//3  # if generation//3 numbers of result are same, then exit the loop
        self.various_name = list(map(lambda x:'x{}'.format(x), list(range(self.various_num))))  # create pandas column name for various dna
        self.main_process()  # main running process
        #___________________________________________________________________

    def main_process(self):
        bar = progressbar.ProgressBar(0, self.generation)
        self.first_chrom()  # generate the first chrom
        self.chrom = self.fchrom.copy()  # shift the first chrom to common various chrom
        self.result = []  # recording result of function for each generation
        for i in range(0, self.generation):
            self.fitness()  # calculate the fitness
            self.result.append(self.chrom.ix[0, 'RESULT'])  # store the best fitness result
            self.elite_storage_fitness()  # storage numbers of elite with best fitness, numbers according to self.elite_fitness_num
            self.relative()  # count the DNA similar condition with other
            self.expectation()  # calculate the expectation
            self.elite_storage_expectation()  # storage the elite expectation
            # store into self.elite_expectation
            self.chrom = self.new_generation()  # create new generation
            bar.update(i)
            # exit loop if all result same within numbers of self.exit_flag
            if i > self.exit_flag:
                if self.result[i] == self.result[i-self.exit_flag]:
                    print()  #print enter
                    print('Result fix in 1/3 of generation, exit loop!')
                    break   # break loop of generation         
        self.fitness()  # calculate the final generation fitness, with RESULT and FITNESS, sorted by fitness
        self.result.append(self.chrom.ix[0, 'RESULT'])  # append the final best result after loop
        best_x = self.mapping(list(self.chrom.ix[0, 0:self.various_num]))  # get the top first dna various and turn to target x range
        plt.plot(list(range(len(self.result))), self.result)  #plot the result in generation
        print()  # print enter
        print('Result as following:\n')
        for n in range(self.various_num):
            print('{x}={y};\n'.format(x=self.various_name[n], y=best_x[n]))
        print('Result={}'.format(self.result[::-1][0]))
        plt.show()
   
    def first_chrom(self):
        for i in range(self.various_num):
            # create first chrom and turn to binary in pandas form: dec to bin, delete '0b', fill zero to DNA length, also create the pandas form
            fchrom = pd.DataFrame([str(bin(random.randint(0, 2**self.dna_length)).replace('0b', '')).zfill(self.dna_length) for i in range(self.population)], columns=[self.various_name[i]])
            if i == 0:
                self.fchrom = fchrom.copy()
            else:
                self.fchrom = pd.concat([self.fchrom, fchrom], axis=1, join='outer')  #join the various dna in one pandas in columns

    def mapping(self, x_value):  # turn dna to dec, and map the DNA range to fitness function range of values
        x_value = list(map(lambda x: int(x, 2), x_value))  # convert bin to dec
        for i in range(self.various_num):
            x_value[i] = (self.upper_edge[i]-self.lower_edge[i])/2**self.dna_length*x_value[i]+self.lower_edge[i]  #  and mapping to the value range of target function
        return x_value

    def mapminmax(self, x):  # mapminnmax, turn list from 0 to 1; x should be a list
        if min(x) == max(x):
            y = list(itertools.repeat(1, len(x)))
        else:
            y = []
            for i in x:
                y.append(1-(i-min(x))/(max(x)-min(x)))  # less result better fitness
        return y

    def fitness(self):  # calculate result and fitness of chrom, and concat to the pandas of chrom
        result = []  # recording result of function
        for index, rows in self.chrom.iterrows(): #  iter all single chrom with various, in row
            x_value = list(rows.ix[0:self.various_num])  # pick up dna of all various and turn to list
            x_value = self.mapping(x_value)  # turn dna to num according to target range
            result.append(self.func(x_value))  # append the result of each various dna to the list 'result'
        fitness = self.mapminmax(result)  # mapminmax result to fitness, max to 1, min to 0
        fitness = pd.DataFrame([result, fitness], index=['RESULT', 'FITNESS'])  #result and fitness to pandas
        fitness = fitness.T  # row to columns, in order to concat
        self.chrom = pd.concat([self.chrom, fitness], axis=1, join='outer')  # concat the result right behind the chrom
        self.chrom = self.chrom.sort_values(["FITNESS"], ascending=False)  # ascending the chrom by fitness from upper to lower
        self.chrom.index = list(range(len(self.chrom.index)))  # reindex after sorting

    def elite_storage_fitness(self):
        # here terms in fitness function already sorted the fitness
        self.elite_fitness = self.chrom.ix[:, 0:self.various_num].head(self.elite_fitness_num).copy()  # store number of top fitness chrom and also the fitenss and result

    def relative(self):
        relative_factor = []  # empty the recording list for all DNA in chrom
        relative_sum_various = []  # to sum relative factor of various in row (for each fitness)
        for i in range(self.various_num):  # to loop all in one various dna, and then next dna
            for dna in self.chrom.ix[:, i]:  # loop, from 1st DNA to last,of each various dna
                relative_num = 1  # similar count recording for one single DNA, value 1 for preventing divition of mapminmax by zero
                for cmp_dna in self.chrom.ix[:, i]:  # loop selected DNA for every other DNA, for counting similar relative, including itself
                    relative_num += le.distance(dna, cmp_dna)  # Levenshtein distance, for counting similar of two strings
                    # one single various dna loop end
                relative_factor.append(relative_num)  # append the counting result for one single DNA of one various
                # population of one various dna loop end
            if relative_sum_various == []:
                relative_sum_various = relative_factor.copy()  # if the first various relative factor, then directively store
            else:
                relative_sum_various = list(map(lambda x, y: x+y, relative_sum_various, relative_factor))  # if not the fist, sum in row for different various in one target
                # all dna loop end
        relative_sum_various = list(map(lambda x: 1 - x/max(relative_sum_various), relative_sum_various))  # mapminmax the relative factor
        # the more similar the value will be more larger (le distance is same to 0, so use 1-x), so to make it as punish factor
        relative_sum_various = pd.DataFrame(relative_sum_various, columns=['RELATIVE'])  # change to pandas
        self.chrom = pd.concat([self.chrom, relative_sum_various], axis=1, join='outer')  # concat to chrom

    def expectation(self):  # calculate the expectation, in order to distinguish a better fitness DNA whether with many relative, if so than lower the possibility
        expectation = self.chrom.FITNESS-self.chrom.RELATIVE*self.relative_valve  # as the script, expectation is pd
        expectation = pd.DataFrame(expectation, columns=['EXPECTATION'])  # turn to pandas,again
        self.chrom = pd.concat([self.chrom, expectation], axis=1, join='outer')  # join to chrom

    def elite_storage_expectation(self): # storage the top numbers of elite of expectation
        self.chrom = self.chrom.sort_values(["EXPECTATION"], ascending=False)  # ascending the chrom by expectation from upper to lower
        self.chrom.index = list(range(len(self.chrom.index)))  # reindex after sorting
        self.elite_expectation = self.chrom.ix[:, 0:self.various_num].head(self.elite_expectation_num).copy()  # store number of top fitness chrom

    def wheel_gambling(self, k):  # wheel gambling, to locate one lucky DNA which prepare to marry :)
        possibility = list(itertools.accumulate(self.chrom.EXPECTATION))  # possibility index of accumulation in expectation
        gamble = random.uniform(0, max(possibility))  # random float number between 0 to max accumulation of possibility
        pos_pin = 0  # a pin to locate the random number drop to which single dna
        for i in possibility:  # search the accumulated expectation
            if gamble <= i:  # when if the random point lower than the accumulate expectation, it means it felt into the single DNA
                break  # then break the searching
            pos_pin += 1  # if not, then pin to next one
        return self.chrom.ix[pos_pin, k]  # return the lucky one, with only DNA columns of the various

    def random_dna_section(self):  #pick up random section of dna, in list of 2 number
        exitflag = True
        while exitflag:  # continue to generate the cross range until they not same
            looks_like = [random.randint(0, self.dna_length+1) for i in range(2)]  # cross range, +1 due to index to include the last one in DNA
            looks_like.sort() #make sure the small one in fore
            if looks_like[0] != looks_like[1]:  # two cross point are not same
                exitflag = False  # then exit loop
        return looks_like

    def cross(self, x):  # cross the lucky dna, dna string in list to x
        sec = self.random_dna_section()  # random section of dna pin, two number in a list
        x_section = []  # for recording dna cross section take out, list for both
        for i in x:
            x_section.append(i[sec[0]:sec[1]])  # dna cross section take out, list for both
        x[0] = x[0].replace(x[0][sec[0]:sec[1]], x_section[1])  # cross, the first dna section replace by section of the second dna
        x[1] = x[1].replace(x[1][sec[0]:sec[1]], x_section[0])  # vice versa
        return x

    def marriage(self, k):  # pickup two lucky dna and generate two kid crossed
        lucky_dna = [self.wheel_gambling(k) for i in range(2)]  # pickup two lucky dna
        baby_dna = self.cross(lucky_dna)  # get two new generation
        return baby_dna

    def new_generation(self):  # create new generation
        loop_generation = int((self.population - self.elite_fitness_num - self.elite_expectation_num)/2)  # time for loop to create new chrom other than the elite
        for k in range(self.various_num):  # due to only one expectation, so loop generation by each various dna
            new_chrom = []  #recorder for new chrom
            for i in range(loop_generation):
                new_generation = self.marriage(k)  # create two cross new chrom
                for j in new_generation:  # assign two to new chrom one by one
                    new_chrom.append(j)
                    # loop for single various dna end
            if k == 0:
                new_chrom_various = pd.DataFrame(new_chrom, columns=[self.various_name[k]])  # new cross chrom for pandas
            else:
                new_chrom = pd.DataFrame(new_chrom, columns=[self.various_name[k]])
                new_chrom_various = pd.concat([new_chrom_various, new_chrom], axis=1, join='outer')
        new_chrom_various = pd.concat([new_chrom_various, self.elite_fitness, self.elite_expectation], axis=0, join='outer')  # concat new chrom, and both elite
        new_chrom_various.index = list(range(len(new_chrom_various.ix[:, 0])))  # reindex the chrom pandas
        new_chrom_various = self.mutation(new_chrom_various)
        return new_chrom_various

    def mutation(self, x):  # mutation, to keep the chrom variable, input x to be chrom pandas
        for i in range(len(x.ix[:, 0])):  # loop all DNA in row
            if random.random() <= self.mutation_possibility: # if random less than the possibility of mutation, than mutation
                for j in range(self.various_num):
                    mutate_dna = []  # to reset the mutate_dna for each loop
                    looks_like = self.random_dna_section()  # mutation section of dna
                    mutate_dna = "".join([random.choice(['0', '1']) for i in range(looks_like[1]-looks_like[0])])  # loop for generate mutate section dna
                    x.ix[i, j] = x.ix[i, j].replace(x.ix[i, j][looks_like[0]:looks_like[1]], mutate_dna)
        return x


myfunc = ga()


