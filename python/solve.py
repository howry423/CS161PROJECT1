"""Solves an instance.

Modify this file to implement your own solvers.

For usage, run `python3 solve.py --help`.
"""

import argparse
from pathlib import Path
from typing import Callable, Dict
from point import Point
import numpy as np
from instance import Instance
from solution import Solution
from file_wrappers import StdinFileWrapper, StdoutFileWrapper
import math
from random import randint, seed


def solve_naive(instance: Instance) -> Solution:
    return Solution(
        instance=instance,
        towers=instance.cities,
    )

def solve_sequential_greedy(instance: Instance) -> Solution:

    # initialize data structures
    solution_set = []

    # set constants
    Rp = instance.R_p
    D = instance.grid_side_length

    # create city matrices and find valid tower positions
    listTower, listCityMatrix, dictCity = create_city_matrices(instance)

    # find highest tower
    while len(dictCity) > 0:
        highestTower = listTower[0]
        connected_cities = []
        highestDegree = 0
        for idx, tower in enumerate(listTower):
            tempDegree = 0
            temp_connected_cities = []
            for city, city_idx in dictCity.items():
                connected = listCityMatrix[city_idx][tower.x][tower.y]
                # if tower is connected to city
                if connected > 0:
                    temp_connected_cities.append(city)
                    tempDegree += connected
            if tempDegree > highestDegree:
                highestDegree = tempDegree
                highestTower = tower
                connected_cities = temp_connected_cities

        solution_set.append(highestTower)
        
        for i in range(highestTower.x - Rp if highestTower.x - Rp >= 0 else 0, highestTower.x + Rp + 1 if highestTower.x + Rp + 1 < D+1 else D):
            for j in range(highestTower.y - Rp if highestTower.y - Rp >= 0 else 0, highestTower.y + Rp + 1 if highestTower.y + Rp + 1 < D+1 else D):
                if Point.distance_obj(Point.parse("{} {}".format(i, j)), highestTower) <= Rp:
                    for cityMatrix in listCityMatrix:
                        cityMatrix[i][j] = cityMatrix[i][j]/(math.exp(0.17))
        # find the cities the tower is connected to and remove it from dictCity
        for city_point in connected_cities:
            dictCity.pop(city_point)

    return Solution(instance=instance, towers=solution_set)

def solve_distributed_greedy(instance: Instance) -> Solution:

    # set seed
    seed(42)

    # initialize data structures
    solution_set = []

    # set constants
    D = instance.grid_side_length
    Rp = instance.R_p

    # create city matrices and find valid tower positions
    listTower, listCityMatrix, dictCity = create_city_matrices(instance)

    # carry out matrix addition to obtain matrix containing degrees of towers
    finalMatrix = sum(listCityMatrix)

    # find highest tower (local maximum)
    while len(dictCity) > 0:

        # select random point from listTower
        rand_index = randint(0, len(listTower)-1)
        highestTower = listTower[rand_index]
        highest_x, highest_y = highestTower.x, highestTower.y
        connected_cities = []
        curDegree = finalMatrix[highest_x][highest_y]

        # check all towers within service radius for highest tower
        while highest_x < D-1 and highest_y < D-1 and highest_x >= 0 and highest_y >= 0:
            for i in range(-2, 3):
                for j in range(-2, 3):
                    x_neighbour = highest_x + i
                    y_neighbour = highest_y + j
                    if (x_neighbour, y_neighbour) in listTower:
                        # check if neighbour tower has higher degree than current tower
                        neighbourDegree = finalMatrix[x_neighbour][y_neighbour]
                        if neighbourDegree > curDegree:
                            highest_x, highest_y = x_neighbour, y_neighbour
                            break
                else:
                    continue
                break        
        
        temp_connected_cities = []
        for city, city_idx in dictCity.items():
            connected = listCityMatrix[city_idx][highest_x][highest_y]
            # if tower is connected to city
            if connected > 0:
                temp_connected_cities.append(city)
        connected_cities = temp_connected_cities

        # add this tower with highest degree to solution set
        solution_set.append((highest_x, highest_y))
        
        # penalizing 
        for i in range(highest_x - Rp if highest_x - Rp >= 0 else 0, highest_x + Rp + 1 if highest_x + Rp + 1 < D+1 else D):
            for j in range(highest_y - Rp if highest_y - Rp >= 0 else 0, highest_y + Rp + 1 if highest_y + Rp + 1 < D+1 else D):
                if Point.distance_obj(Point.parse("{} {}".format(i, j)), (highest_x, highest_y)) <= Rp:
                    for cityMatrix in listCityMatrix:
                        cityMatrix[i][j] = cityMatrix[i][j]/(math.exp(0.17))

        # find the cities the tower is connected to and remove it from dictCity
        for city_point in connected_cities:
            dictCity.pop(city_point)

    return Solution(instance=instance, towers=solution_set)

# HELPER FUNCTION
def create_city_matrices(instance: Instance):
    
    # initialize data structures
    listCityMatrix = []
    dictCity = {}
    listTower = []
    setTower = []

    # set constants
    D = instance.grid_side_length
    Rp = instance.R_p

    # for loop to fill up matrix for each city
    for idx, city in enumerate(instance.cities):
        dictCity[city] = idx
        cityMatrix = np.zeros((D, D))
        xCoord = city.x
        yCoord = city.y
        for i in range(-2, 3):
            for j in range(-2, 3):
                xresult = xCoord + i
                yresult = yCoord + j
                xresult = xresult if xresult >= 0 else 0
                xresult = xresult if xresult < D else D-1
                yresult = yresult if yresult >= 0 else 0
                yresult = yresult if yresult < D else D-1
                cityMatrix[xresult][yresult] = 1
                setTower.append((xresult, yresult))
        xresult = xCoord + 3 if xCoord + 3 < D else D-1
        cityMatrix[xresult][yCoord] = 1
        setTower.append((xresult, yresult))
        xresult = xCoord - 3 if xCoord - 3 >= 0 else 0
        cityMatrix[xresult][yCoord] = 1
        setTower.append((xresult, yresult))
        yresult = yCoord + 3 if yCoord + 3 < D else D-1
        cityMatrix[xCoord][yresult] = 1
        setTower.append((xresult, yresult))
        yresult = yCoord - 3 if yCoord - 3 >= 0 else 0
        cityMatrix[xCoord][yresult] = 1
        setTower.append((xresult, yresult))
        listCityMatrix.append(cityMatrix)
    setTower = list(set(setTower))
    for tower in setTower:
        listTower.append(Point.parse("{} {}".format(tower[0], tower[1])))
    
    return listTower, listCityMatrix, dictCity


SOLVERS: Dict[str, Callable[[Instance], Solution]] = {
    "naive": solve_naive,
    "greedy": solve_sequential_greedy,
    "greedy_distributed": solve_distributed_greedy
}


# You shouldn't need to modify anything below this line.
def infile(args):
    if args.input == "-":
        return StdinFileWrapper()

    return Path(args.input).open("r")


def outfile(args):
    if args.output == "-":
        return StdoutFileWrapper()

    return Path(args.output).open("w")


def main(args):
    with infile(args) as f:
        instance = Instance.parse(f.readlines())
        solver = SOLVERS[args.solver]
        solution = solver(instance)
        assert solution.valid()
        with outfile(args) as g:
            print("# Penalty: ", solution.penalty(), file=g)
            solution.serialize(g)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve a problem instance.")
    parser.add_argument("input", type=str, help="The input instance file to "
                        "read an instance from. Use - for stdin.")
    parser.add_argument("--solver", required=True, type=str,
                        help="The solver type.", choices=SOLVERS.keys())
    parser.add_argument("output", type=str,
                        help="The output file. Use - for stdout.",
                        default="-")
    main(parser.parse_args())
