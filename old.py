

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
    network_directory: str, network_file: str = "original_network.lmp"
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
    original_network_file = os.path.realpath(
        os.path.join(network_directory, network_file)
    )
    network_file = os.path.realpath(os.path.join(network_directory, "network.lmp"))
    output_network_file = os.path.realpath(
        os.path.join(network_directory, "output_network.lmp")
    )
    uncompressed_network = Network.from_data_file(
        original_network_file, include_angles=False, include_dihedrals=False
    )
    uncompressed_network.write_to_file(network_file)
    run_lammps(network_directory, mode="single", num_threads=1, num_procs=1)
    compressed_network = Network.from_data_file(
        output_network_file, include_angles=False, include_dihedrals=False
    )
    transverse_strain = compressed_network.box.x - uncompressed_network.box.x
    axial_strain = compressed_network.box.y - uncompressed_network.box.y
    return p_ratio(transverse_strain, axial_strain)
