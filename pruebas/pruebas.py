#!/usr/bin/env python3

"""
    ModelBark: a toy model to study bark formation in Woody species

    Authors: Álvaro Gutiérrez Climent, Juan Carlos Nuño, Unai López de Heredia, Álvaro Soto

    UTF-8

"""


import tkinter as tk
import os
from PIL import ImageTk, Image
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

working_directory = os.getcwd()

class Plant:
    def __init__(self):

        self.radius = [1]
        self.number_phellogen = 0
        self.xylem_storage = [0]
        self.phloem_storage = [0]
        self.inactive_phloem_storage = [0]
        self.phellem_storage = [0]
        self.phelloderm_storage = [0]
        self.equation_storage = [1]
        self.first_phellogen_storage = 0

    def sample(self, probability: float):

        sample = np.random.random_sample()

        if sample <= probability:
            return int(0)

        else:
            return int(1)

    def vascular_cambium_division(self, vcdr: float):

        radius_mod = self.radius

        sample = self.sample(vcdr)

        if sample == 0:
            radius_mod.insert(0, 0)

        if sample == 1:
            radius_mod.insert(radius_mod.index(1) + 1, 0)

        self.radius = radius_mod

    def first_phellogen(self, percentage: float):

        radius_mod = self.radius

        vascular_cambium_position = radius_mod.index(1)

        radius_length = len(radius_mod)

        insert_place = (
            radius_length - (round((radius_length - vascular_cambium_position) * percentage)))

        radius_mod.insert(insert_place, 3)

        self.radius = radius_mod

        self.first_phellogen_storage = len(radius_mod)

        self.number_phellogen = 1

    def phellem_production(self, pdr: float):

        radius_mod = self.radius

        number_phellogen_mod = self.number_phellogen

        phellogen_position = radius_mod.index(3)

        if self.sample(pdr) == 0:
            radius_mod.insert(phellogen_position + 1,
                              (number_phellogen_mod + 3))

        self.radius = radius_mod

    def phelloderm_production(self, pddr: float):

        radius_mod = self.radius

        phellogen_position = radius_mod.index(3)

        if self.sample(pddr) == 0:
            radius_mod.insert(phellogen_position,
                              (2))
        
        self.radius = radius_mod

        
    def new_phellogen(self, percentage: float):

        radius_mod = self.radius

        vascular_cambium_position = radius_mod.index(1)

        phellogen_position = radius_mod.index(3)

        insert_place = (phellogen_position - round((phellogen_position -
                        vascular_cambium_position) * percentage))

        radius_mod[phellogen_position] = (self.number_phellogen + 3)

        radius_mod.insert(insert_place, 3)

        self.number_phellogen += 1

        self.radius = radius_mod

    def radius_length(self):

        return len(self.radius)

    def num_last_phellogen_cells(self):

        radius_mod = self.radius

        number_phellogen = self.number_phellogen

        return int(radius_mod.count(number_phellogen + 3))

    def num_xylem_cells(self):

        radius_mod = self.radius

        return radius_mod.index(1)

    def num_phloem_cells(self):

        radius_mod = self.radius

        vascular_cambium_position = radius_mod.index(1)

        if radius_mod.count(3) == 1:

            phellogen_position = radius_mod.index(3)

            return radius_mod[vascular_cambium_position:phellogen_position].count(0)
            
        else:
            return radius_mod.count(0) - vascular_cambium_position

    def num_phellem_cells(self):

        radius_mod = self.radius

        if radius_mod.count(3) == 0:
            return 0

        else:

            num_xylem_or_phloem_cells = radius_mod.count(0)
            num_phelloderm_cells = radius_mod.count(2)
            radius_length = len(radius_mod)
            vascular_suber_cambium_cells = 2
            return radius_length - num_xylem_or_phloem_cells - num_phelloderm_cells - vascular_suber_cambium_cells
        
        
    def num_phelloderm_cells(self):

        radius_mod = self.radius

        return radius_mod.count(2)


    def num_inactive_phloem_cells(self):

        radius_mod = self.radius

        num_xylem_or_phloem_cells = radius_mod.count(0)

        xylem_cells = self.num_xylem_cells()

        active_phloem_cells = self.num_phloem_cells()

        return num_xylem_or_phloem_cells - xylem_cells - active_phloem_cells

    def parameters(self):

        xylem = self.num_xylem_cells()

        phloem = self.num_phloem_cells()

        phellem = self.num_phellem_cells()

        inactive_phloem = self.num_inactive_phloem_cells()

        phelloderm = self.num_phelloderm_cells()

        return [xylem, phloem, phellem, inactive_phloem, phelloderm]

    def equation(self, a, b, c, d, e):

        radius_parameters = self.parameters()

        xylem_a = a * radius_parameters[0]

        phloem_b = b * radius_parameters[1]

        phellem_c = c * radius_parameters[2]

        inactive_phloem_d = d * radius_parameters[3]

        phelloderm_e = e * radius_parameters[4]

        return (1 + phellem_c + inactive_phloem_d + phelloderm_e) / (xylem_a + phloem_b)

    def graphical_parameters_storage(self, a, b, c, d, e):

        self.xylem_storage.append(self.num_xylem_cells())

        self.phloem_storage.append(self.num_phloem_cells())

        self.phellem_storage.append(self.num_phellem_cells())

        self.inactive_phloem_storage.append(self.num_inactive_phloem_cells())

        self.phelloderm_storage.append(self.num_phelloderm_cells)

        self.equation_storage.append(self.equation(a, b, c, d, e))

    def result(self):

        output = self.parameters()

        output.append(self.number_phellogen)

        output.append(self.first_phellogen_storage)

        output.append(self.radius)

        return output

    def plotting(self, k):

        if not os.path.exists(os.path.join(working_directory, 'Figures')):
            os.mkdir(os.path.join(working_directory, 'Figures'))

        # Parameters through time plot

        plt.plot(self.xylem_storage, label='Xylem')

        plt.plot(self.phloem_storage, label='Phloem')

        plt.plot(self.phellem_storage, label='Phellem')

        plt.plot(self.inactive_phloem_storage, label='Inactive Phloem')

        plt.plot(self.phelloderm_storage, label= 'Phelloderm')

        plt.legend()

        plt.title('Parameters of the model through time')

        plt.savefig('Figures/Parameters_plot.jpg')

        plt.close()

        # Equation through time plot

        plt.plot(self.equation_storage, label='F(t)')

        plt.axhline(y=k, color='r', linestyle=':', label='K threshold')

        plt.legend()

        plt.title('F(t) through time')

        plt.savefig('Figures/Equation_plot.jpg')

        plt.close()


def simulation_generation(vascular_cambium_division_rate: float, phellogen_division_rate: float, phelloderm_division_rate: float ,phellogen_position: float, a: float, b: float, c: float, d: float, e:float, threshold: float, max_length: int):
    simulation = Plant()
    simulation.vascular_cambium_division(vascular_cambium_division_rate)
    simulation.graphical_parameters_storage(a, b, c, d, e)

    while simulation.equation(a, b, c, d, e) >= threshold:
        simulation.vascular_cambium_division(vascular_cambium_division_rate)
        simulation.graphical_parameters_storage(a, b, c, d, e)

        if simulation.radius_length() >= max_length:
            break

    if simulation.radius_length() >= max_length:
        simulation.plotting(threshold)
        return simulation.result()

    else:

        simulation.first_phellogen(phellogen_position)
        while simulation.radius_length() <= max_length:

            while simulation.equation(a, b, c, d, e) >= threshold and simulation.radius_length() <= max_length:
                simulation.vascular_cambium_division(
                    vascular_cambium_division_rate)
                simulation.phellem_production(phellogen_division_rate)
                simulation.phelloderm_production(phelloderm_division_rate)
                simulation.graphical_parameters_storage(a, b, c, d, e)

            if simulation.radius_length() <= max_length:
                simulation.new_phellogen(phellogen_position)
                simulation.vascular_cambium_division(
                    vascular_cambium_division_rate)

            else:
                break

        # simulation.plotting(threshold)
        return simulation.result()
    

print(simulation_generation(0.9,0.01,0.1,0.3,0.016,0.008,0.8,0.002,0.001,1,1000))



