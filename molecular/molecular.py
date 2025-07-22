from qiskit_nature.drivers import Molecule
from qiskit_nature.drivers.second_quantization import ElectronicStructureMoleculeDriver, ElectronicStructureDriverType
from qiskit_nature.problems.second_quantization import ElectronicStructureProblem
from qiskit_nature.transformers.second_quantization.electronic import ActiveSpaceTransformer
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.mappers.second_quantization import JordanWignerMapper


def LiH(dist):
    molecule = Molecule(
        geometry=[["li", [0, 0, .0]], ["H", [dist, 0, .0]]],
        charge=0, multiplicity=1,
    )
    driver = ElectronicStructureMoleculeDriver(
        molecule, basis="6-31g", driver_type=ElectronicStructureDriverType.PYSCF,
    )
    es_problem = ElectronicStructureProblem(driver,[ActiveSpaceTransformer(2,6)])
    second_q_op = es_problem.second_q_ops()
    qubit_converter = QubitConverter(mapper=JordanWignerMapper())
    qubitOp = qubit_converter.convert(second_q_op['ElectronicEnergy'])
    H=qubitOp.to_spmatrix()
    return H,qubitOp 