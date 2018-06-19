import pyrosetta
from pyrosetta import rosetta


def get_task_factory(pose, designable_residues, repackable_residues, extra_rotamers=True, limit_aro_chi2=True, layered_design=True):
    '''Get a task factory given the designable and repackable residues.'''
    def list_to_str(l):
        return ','.join(list(str(i) for i in l))

    task_factory = rosetta.core.pack.task.TaskFactory()

    if len(designable_residues) > 0:
        designable_selector = rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str(designable_residues)) 
        racaa = rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        racaa.aas_to_keep('GAPVILMFYWSTKRDENQ') # No CYS or HIS
        designable_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                racaa, designable_selector)
        task_factory.push_back(designable_operation)

    if len(repackable_residues) > 0:
        repackable_selector = rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str(repackable_residues)) 
        repackable_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                rosetta.core.pack.task.operation.RestrictToRepackingRLT(), repackable_selector)
        task_factory.push_back(repackable_operation)

    natro_residues = [i for i in range(1, pose.size() + 1) if not i in designable_residues + repackable_residues]
    if len(natro_residues) > 0:
        natro_selector = rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str(natro_residues)) 
        natro_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                rosetta.core.pack.task.operation.PreventRepackingRLT(), natro_selector)
        task_factory.push_back(natro_operation)

    if extra_rotamers:
        ers = rosetta.core.pack.task.operation.ExtraRotamersGeneric()
        ers.ex1(True)
        ers.ex2(True)
        ers.extrachi_cutoff(0)

        task_factory.push_back(ers)

    if limit_aro_chi2:
        lac = rosetta.protocols.task_operations.LimitAromaChi2Operation()

        task_factory.push_back(lac)

    if layered_design:
        ld = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_task_operation(
            '''<LayerDesign name="layer_all" layer="core_boundary_surface_Nterm_Cterm" use_sidechain_neighbors="True">
    		<Nterm>
    			<all append="DEGHKNQRST" />
    			<all exclude="CAFILMPVWY" />
    		</Nterm>
    		<Cterm>
    			<all append="DEGHKNQRST" />
    			<all exclude="CAFILMPVWY" />
    		</Cterm>
        </LayerDesign>''')
        task_factory.push_back(ld)

    return task_factory

