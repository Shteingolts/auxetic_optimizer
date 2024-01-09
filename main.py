from __future__ import annotations

import dataclasses
import os
import pickle
import shutil
import subprocess
import sys
import time
from copy import deepcopy
from multiprocessing import Pool

from network import Bond, Network


def check_lammps_installation() -> str | None:
    """Checks wheter lammps is installed on the current machine

    Returns
    -------
    str | None
        Either path to the lammps executable or None if lammps is not found
    """
    try:
        out = subprocess.check_output(["which", "lmp"]).decode("utf-8")
        print(f"Lammps found: {out}")
        return out
    except subprocess.CalledProcessError:
        print("Lammps is probably not installed on your machine")
        return None


def clean_up(
    calculation_directory: str, network_filename: str = "original_network.lmp"
):
    """
    A helper function which deletes all the angle and dihedral data from the file.
    """
    network_filepath = os.path.join(calculation_directory, network_filename)
    network = Network.from_data_file(
        network_filepath, include_angles=False, include_dihedrals=False
    )
    network.write_to_file(network_filepath)


def p_ratio(transverse_strain: float, axial_strain: float) -> float:
    """Poisson's ratio - measure of Poisson's effect, i.e.,
    deformation in the direction perpendicular to direction
    of loading. Elongating dimension devided by shrinking dimension.

    Parameters
    ----------
    trans_strain :  
        Increasing dimension
    axial_strain : float
        Shrinking dimension

    Returns
    -------
    float
        P ratio
    """
    return -axial_strain / transverse_strain


def p_ratio_from_file(
        network_directory: str,
        network_file: str = "original_network.lmp"
) -> float:
    """Computes Poisson's ratio from a network file by running compression simulation.

    Parameters
    ----------
    network_directory : str
        Directory containing a network file and simulation files (see example folder)
    network_file : str, optional
        filename for the network file, by default "original_network.lmp"

    Returns
    -------
    float
        poisson's ratio
    """
    original_network_file = os.path.realpath(os.path.join(network_directory, network_file))
    network_file = os.path.realpath(os.path.join(network_directory, "network.lmp"))
    output_network_file = os.path.realpath(os.path.join(network_directory, "output_network.lmp"))
    uncompressed_network = Network.from_data_file(original_network_file, include_angles=False, include_dihedrals=False)
    uncompressed_network.write_to_file(network_file)
    run_lammps(network_directory, mode='single', num_threads=1, num_procs=1)
    compressed_network = Network.from_data_file(output_network_file, include_angles=False, include_dihedrals=False)
    transverse_strain = compressed_network.box.x - uncompressed_network.box.x
    axial_strain = compressed_network.box.y - uncompressed_network.box.y
    return p_ratio(transverse_strain, axial_strain)    


def run_lammps(
    calculation_directory: str,
    mode: str = "single",
    num_threads: int = 1,
    num_procs: int = 1,
):
    """
    A helper function which runs the external lammps code.
    """
    input_file_path = os.path.join(calculation_directory, "in.compression")
    mpi_command = f"mpirun -np {num_procs} lmp -in in.stress".split()
    command = f"lmp -in {input_file_path}".split()

    os.chdir(os.path.dirname(input_file_path))
    os.environ["OMP_NUM_THREADS"] = str(num_threads)
    if mode == "single":
        subprocess.run(command, stdout=subprocess.DEVNULL)
    elif mode == "mpi":
        subprocess.run(mpi_command, stdout=subprocess.DEVNULL)


@dataclasses.dataclass
class CalculationResult:
    """Dataclass which combines the result of the network optimization step.
    """
    bond: Bond
    p_ratio: float
    others: list[CalculationResult] | None = None

    def __repr__(self) -> str:
        return f"Calc result: Bond: {self.bond}, P ratio = {self.p_ratio}"


@dataclasses.dataclass
class CalculationSetup:
    """
    A dataclass accepted by `run_iteration()` function.
    A collection of arguments basically.
    """
    network_directory: str
    bonds_list: list[Bond]
    network_file: str = "original_network.lmp"


class StepResult():
    step_number: int
    p_ratio: float
    bond_removed: Bond
    network: Network
    closest_contenders: None | list[CalculationResult]

    def __init__(
            self,
            step_number: int,
            p_ratio: float,
            bond_removed: Bond,
            network: Network,
            closest_contenders: None | list[CalculationResult] = None
    ) -> None:
        self.step_number = step_number
        self.network = network
        self.p_ratio = p_ratio
        self.bond_removed = bond_removed
        self.closest_contenders = closest_contenders
    
    def __repr__(self) -> str:
        return f"Step {self.step_number}: FF={self.p_ratio}, Bond removed: {self.bond_removed}"


def run_iteration(setup: CalculationSetup) -> CalculationResult:
    """
    Deletes every bond from the network once, then runs lammps compression simulation
    to calculate the Poisson's ratio. Returns the best bond to delete and all the others that had a
    positive impact, i.e., decreased the P ratio.
    Needs to have a file named `original_network.lmp` in the provided directory.

    Parameters
    ----------
    network_directory : str
        Directory containing a network file and simulation files (see example folder)
    network_file : str, optional
        filename for the network file, by default "original_network.lmp"

    Returns
    -------
    CalculationResult
        
    """

    initial_p = p_ratio_from_file(setup.network_directory)
    original_network_file = os.path.join(setup.network_directory, setup.network_file)
    output_network_file = os.path.join(setup.network_directory, "output_network.lmp")
    print(setup.network_directory)
    print(original_network_file)
    print(output_network_file)

    original_network = Network.from_data_file(original_network_file, include_angles=False, include_dihedrals=False)

    current_p: float = initial_p
    bond_to_be_deleted: Bond | None = None
    others: list[Bond] = []
    for index, bond in enumerate(setup.bonds_list):
        
        # deepcopy because we don't want to modify the original network
        current_network = deepcopy(original_network)

        # remove a bond and write the network to a file
        current_network.remove_bond(bond)
        current_network.write_to_file("network.lmp")

        # run lammps on the previously written network
        run_lammps(setup.network_directory, mode="single")

        # read the minimized network
        # calculate source strain, target strains, and fitness function
        compressed_network = Network.from_data_file(output_network_file, include_angles=False, include_dihedrals=False)
        transverse_strain = compressed_network.box.x - current_network.box.x
        axial_strain = compressed_network.box.y - current_network.box.y
        new_p = p_ratio(transverse_strain, axial_strain)

        if new_p < current_p and len(bond.atom1.bonded) >= 4 and len(bond.atom2.bonded) >= 4:
            current_p = new_p
            print(f"New P! {new_p}")
            bond_to_be_deleted = bond
        elif new_p < initial_p and len(bond.atom1.bonded) >= 4 and len(bond.atom2.bonded) >= 4:
            others.append(CalculationResult(bond, new_p, None))
            print(f"Meh P {new_p}")
        else:
            print(f"Worse p {new_p}")
            pass
            # don't need to do anything when the removal 
            # of the bond results in the increase of P ratio
    return CalculationResult(bond_to_be_deleted, current_p, others)


def load_optimization_log(pickle_file: str) -> list[StepResult]:
    optimization_log: list[StepResult] = []
    with open(pickle_file, 'rb') as log_file:
        try:
            while True:
                optimization_log.append(pickle.load(log_file))
        except EOFError:
            pass
    return optimization_log


def print_history(log_file_path: str, n_closest: int = 5):
    history = load_optimization_log(log_file_path)
    for step in history[:10]:
        print(f"{step.step_number}: {round(step.p_ratio, 4)}")
        for index, close in enumerate(step.closest_contenders[:5]):
            print(f"  +{round(close.p_ratio-step.p_ratio, 4)}")
        print("--------------------")


def parallel(
    main_calculation_directory: str,
    n_procs: int = 6,
    p_threshold: float = 0.01,
) -> dict[int, StepResult]:
    """
    The main function which performs the network optimization.

    For every iterative step of the optimization, a new directory named "step_N" is created,
    where N is the serial number for the current iteration, starting from 1.
    Multiprocessing is implemented in way such that for every core requested a seperate directory
    named "core_N" is created within the "step_N" directory.
    The total number of bond deletions are separated into M batches, where M=n_bonds/_ncores. Each core 
    is assigned its own batch. After all indivial iteration are complete, one bond is removed from the network which resulted
    in the highest decrease of fitness function.

    The process runs until either:
    - 1. There is no change in fitness function observed between iterations, or
    - 2. Fitness function reaches a certain threshold difined by `ff_threshold` variable.

    The number of core available for the process is controled by the `n_procs` variable.
    """
    start_time = time.perf_counter()
    step_counter = 1
    intermidiate_results = {}

    # create a log file. Overwrite if it's there already
    log_file_path = os.path.join(main_calculation_directory, 'optimization_log.pkl')
    log_file = open(log_file_path, 'wb')
    log_file.close()
    while True:
        # create a directory for the current calculation step
        current_working_directory = os.path.realpath(os.path.join(main_calculation_directory, f"step_{step_counter}"))
        print(f"cwd {current_working_directory}")
        os.makedirs(current_working_directory)

        # if it's not the first iteration, copy `result.lmp` file 
        # from the previous step directory and rename it accordingly
        print(f"Step {step_counter}")
        if step_counter != 1:
            shutil.copy(
                os.path.join(main_calculation_directory, f"step_{step_counter-1}", "result.lmp"),
                os.path.join(current_working_directory)
            )
            os.rename(
                os.path.join(current_working_directory, "result.lmp"),
                os.path.join(current_working_directory, "original_network.lmp"),
            )

        # divide all the bonds in the network into a number of tasks equal to the number of cores
        network = (Network.from_data_file(os.path.join(main_calculation_directory, "original_network.lmp"), include_angles=False, include_dihedrals=False)
            if step_counter == 1
            else Network.from_data_file(
                os.path.join(current_working_directory, "original_network.lmp")
            ))
        
        n_tasks = len(network.bonds) // n_procs
        remainder = len(network.bonds) % n_procs
        tasks = []
        start = 0
        for i in range(n_procs):
            size = n_tasks + (i < remainder)
            tasks.append(network.bonds[start : start + size])
            start += size

        # for each task, create a separate directory within the step directory
        # copy the calculation files
        # create a list of calculation setups
        calc_files = [
            os.path.join(main_calculation_directory, "in.compression"),
            os.path.join(main_calculation_directory, "init.mod"),
            os.path.join(main_calculation_directory, "potential.mod"),
        ]

        calculations: list[CalculationSetup] = []
        for i, task in enumerate(tasks):
            core_dir = os.path.realpath(os.path.join(current_working_directory, f"core_{i+1}"))
            os.makedirs(core_dir)
            # print(core_dir)
            calculations.append(CalculationSetup(core_dir, task))
            for file in calc_files:
                shutil.copy(file, core_dir)

        # copy the network file into a new directory
        if step_counter == 1:
            for i in range(len(tasks)):
                shutil.copy(
                    os.path.join(main_calculation_directory, "original_network.lmp"),
                    os.path.join(current_working_directory, f"core_{i+1}"),
                )
                shutil.copy(
                    os.path.join(main_calculation_directory, "original_network.lmp"),
                    os.path.join(current_working_directory),
                )
        else:
            # also copy the same file into every directory
            for i in range(len(tasks)):
                shutil.copy(
                    os.path.join(main_calculation_directory, f"step_{step_counter-1}", "result.lmp"),
                    os.path.join(current_working_directory, f"core_{i+1}"),
                )
                os.rename(
                    os.path.join(
                        current_working_directory, f"core_{i+1}", "result.lmp"
                    ),
                    os.path.join(
                        current_working_directory, f"core_{i+1}", "original_network.lmp"
                    ),
                )

        with Pool(n_procs) as pool:
            results: list[CalculationResult] = pool.map(run_iteration, calculations)

        # filter out results where the deletion of a bond did not improve the fitness function
        results = [result for result in results if result.bond is not None]

        if results:
            # sort list of CalculationResults by the value of fitness function from lowest to highest    
            results = sorted(results, key=lambda result: result.p_ratio)
            best_result = results[0]
            others: list[CalculationResult] = []
            for index, result in enumerate(results):
                if index == 0 and len(result.others) > 0:
                    for other in result.others:
                        others.append(other)
                else:
                    others.append(CalculationResult(result.bond, result.p_ratio))
                    for other in result.others:
                        others.append(other)
            
            others = sorted(others, key=lambda result: result.p_ratio)
            print(f"Step {step_counter} => {format(best_result.p_ratio, '.3f')}")

            if step_counter == 1:
                okay_to_continue = True
            else:
                okay_to_continue = intermidiate_results[step_counter-1].p_ratio - best_result.p_ratio >= p_threshold

            if not okay_to_continue:
                end_time = time.perf_counter()
                print(f"Threshold {p_threshold} achieved at step {step_counter}.\nElapsed time: {format(start_time - end_time, '.3f')} s")
                network = Network.from_data_file(os.path.join(current_working_directory, "original_network.lmp"))
                network.remove_bond(best_result.bond)
                intermidiate_results[step_counter] = StepResult(step_counter, best_result.p_ratio, best_result.bond, network, closest_contenders=others)
                final_network_file = os.path.join(main_calculation_directory, "final_result.lmp")
                network.write_to_file(final_network_file)
                print(f"The optimized network was written in {final_network_file}")

                with open(log_file_path, 'ab') as log_file:
                    pickle.dump(StepResult(step_counter, best_result.p_ratio, best_result.bond, network, closest_contenders=others), log_file)
                break
            else:
                network = Network.from_data_file(os.path.join(current_working_directory, "original_network.lmp"))
                network.remove_bond(best_result.bond)
                intermidiate_results[step_counter] = StepResult(step_counter, best_result.p_ratio, best_result.bond, network, closest_contenders=others)
                network.write_to_file(os.path.join(current_working_directory, "result.lmp"))
                with open(log_file_path, 'ab') as log_file:
                    pickle.dump(StepResult(step_counter, best_result.p_ratio, best_result.bond, network, closest_contenders=others), log_file)
                step_counter += 1
        else:  
            # if no improvement was observed at all, end the process
            end_time = time.perf_counter()
            print(f"No change in fitness function at step {step_counter}.\nElapsed time: {format(start_time - end_time, '.3f')} s")
            break

    return intermidiate_results


# def main():
#     # usage_info = "\n[USAGE]:\n\n    python3 ./optimize.py\n"
    
#     print("Calculation directory: ")
#     print(">", end='')
#     calculation_directory = input()
#     print("Network file: ")
#     print(">", end='')
#     network_file = input()
#     if os.path.isdir(calculation_directory) and len(os.listdir(calculation_directory)) >= 4:
#         print(f"Calculation directory: {os.path.abspath(calculation_directory)}")
#         calculation_files = {
#             "in.stress" : 0,
#             "init.mod": 0,
#             "potential.mod" : 0,
#         }
#         for file in os.listdir(calculation_directory):
#             if file.strip() in calculation_files:
#                 calculation_files[file.strip()] = 1
        
#         for file, count in calculation_files.items():
#             if count > 0:
#                 print(f"File {file} was found")
#             else:
#                 print(f"File {file} is not found!")
#                 sys.exit(0)
#     else:
#         print(f"{calculation_directory} probably not a valid path")
#         sys.exit(0)

#     print("Source beads (2 integers separated by space): ")
#     print(">", end='')
#     source_beads = input()
#     source_beads = source_beads.strip().split()
#     try:
#         first = int(source_beads[0])
#         second = int(source_beads[1])
#         original_source_beads = (first, second)
#     except ValueError or IndexError:
#         print('Something went wrong.')
#         sys.exit(0)

#     print("Target beads (2 integers separated by space): ")
#     print(">", end='')
#     target_beads = input()
#     target_beads = target_beads.strip().split()
#     try:
#         first = int(target_beads[0])
#         second = int(target_beads[1])
#         original_target_beads = (first, second)
#     except ValueError or IndexError:
#         print('Something went wrong.')
#         sys.exit(0)

#     print("Number of cores: ")
#     print(">", end='')        
#     n_cores = input()
#     try:
#         n_cores = int(n_cores)
#     except ValueError:
#         print(f"{n_cores} is not a digid.")
#         sys.exit(0)
    
#     input_network = os.path.join(calculation_directory, network_file)
#     source_beads, target_beads = prepare_network(
#         input_network,
#         calculation_directory,
#         original_source_beads,
#         original_target_beads
#     )
#     parallel(calculation_directory, source_beads, target_beads, n_procs=n_cores)


if __name__ == '__main__':
    network_directory = os.path.realpath("example")
    parallel(network_directory, n_procs=4, p_threshold=0.002)