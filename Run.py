
import os
from anytree import RenderTree
from vtk.util.numpy_support import vtk_to_numpy
import vtk
from Read_XPLT_Funcs import GetFEB, xpltObj, GetMeshInfo, GetData
import numpy as np

#Choose xplt file
filename = 'resample.xplt'

'''
Construct Tree and get key info
    feb - data tree
    file_size - data size of file
    nstates - nnumber of states
'''
feb,file_size,nstates = GetFEB(filename)

#Draw tree sructure in terminal
for pre,fill,node in RenderTree(feb):
    treestr = u"%s%s" %(pre,node.name)
    print(treestr.ljust(8),node.name,node.read)

'''
Construct Tree and get key info
    nNodes - number of nodes
    nElems - number of elements
    nVar   - nnumber of variables
    StateTimes - time of each state
    VarNames - The names of the variables (there names will be used as the array names later)
    VarType - the types of the variables (respective to the VarNames)
                0 = single precision (s.p.) floating point
                1 = 3D vector of s.p. floats (stored in x, y, z order)
                2 = symmetric 2nd order tensor of s.p. floats (stored in xx, yy, zz, xy, yz, xz order)
                (See FeBio documentation for further info)
    '''
nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
print('The Variable Names are:')
for i in range(nVar):
    print(VarNames[i])

#Collect the data, with all states, note spaces in vairable names are replaced with _
for i in range(nVar):
    VN = VarNames[i].replace(" ","_")
    exec(VN+' = GetData(feb,VarNames[i],nstates,nVar)')
  
#Gather stresses and displacements Note this is dependent on the Var Type
Stress_X  = np.zeros((nstates,nElems))
Stress_Y  = np.zeros((nstates,nElems))
Stress_Z  = np.zeros((nstates,nElems))
Stress_XY = np.zeros((nstates,nElems))
Stress_YZ = np.zeros((nstates,nElems))
Stress_XZ = np.zeros((nstates,nElems))
Disp_X    = np.zeros((nstates,nNodes))
Disp_Y    = np.zeros((nstates,nNodes))
Disp_Z    = np.zeros((nstates,nNodes))
Disp_FEB  = np.zeros((nstates,nNodes,3))      
for i in range(nstates): 
    for j in range(nNodes):
        Disp_X[i,j]  = displacement[i][j*3+1]
        Disp_Y[i,j]  = displacement[i][j*3+2]
        Disp_Z[i,j]  = displacement[i][j*3+3]
    Disp_FEB[i,:,0]  = Disp_X[i,:]  
    Disp_FEB[i,:,1]  = Disp_Y[i,:]  
    Disp_FEB[i,:,2]  = Disp_Z[i,:]
    
for i in range(nstates): 
    for j in range(nElems):
        Stress_X[i,j]  = stress[i][j*6+1]
        Stress_Y[i,j]  = stress[i][j*6+2]
        Stress_Z[i,j]  = stress[i][j*6+3]
        Stress_XY[i,j] = stress[i][j*6+4]
        Stress_YZ[i,j] = stress[i][j*6+5]
        Stress_XZ[i,j] = stress[i][j*6+6]
        
        
reader = vtk.vtkPolyDataReader()

#Save a new vtk file for each state
for j in range(nstates):
    # Read the source file.
    reader.SetFileName('resample.vtk')
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    polydata = reader.GetOutput()
    
    # Add Cell Data on Points
    CellData = [Stress_X[j], Stress_Y[j], Stress_Z[j], Stress_XY[j], Stress_YZ[j], Stress_XZ[j]]
    CellNames = ['Stress_X','Stress_Y','Stress_Z','Stress_XY','Stress_YZ','Stress_XZ']
    
    for i in range(len(CellNames)) :
        arrayCell = vtk.util.numpy_support.numpy_to_vtk(CellData[i], deep=True)
        arrayCell.SetName(CellNames[i])
        dataCells = polydata.GetCellData()
        dataCells.AddArray(arrayCell)
        dataCells.Modified()
        
    # Convert Cell Data to Point Data
    c2p = vtk.vtkCellDataToPointData()
    c2p.AddInputData(polydata)
    c2p.Update()
    c2p.GetOutput()
    NumOfArr = c2p.GetPolyDataOutput().GetPointData().GetNumberOfArrays()
    
    for i in range(NumOfArr):
        if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_X':
            Stress_X_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
        if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_Y':
            Stress_Y_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
        if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_Z':
            Stress_Z_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
        if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_XY':
            Stress_XY_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
        if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_YZ':
            Stress_YZ_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
        if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_XZ':
            Stress_XZ_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
    
    # Add Point Data
    PointData = [Stress_X_Pt,Stress_Y_Pt,Stress_Z_Pt,Stress_XY_Pt,Stress_YZ_Pt,Stress_XZ_Pt,Disp_X[j],Disp_Y[j],Disp_Z[j]]
    PointNames = ['Stress_X_Pt','Stress_Y_Pt','Stress_Z_Pt','Stress_XY_Pt','Stress_YZ_Pt','Stress_XZ_Pt','Disp_X','Disp_Y','Disp_Z']
    
    for i in range(len(PointNames)) :
        arrayPoint = vtk.util.numpy_support.numpy_to_vtk(PointData[i], deep=True)
        arrayPoint.SetName(PointNames[i])
        dataPoints = polydata.GetPointData()
        dataPoints.AddArray(arrayPoint)
        dataPoints.Modified()
    
    # Add Vector Data on Points
    VectorData = [Disp_FEB[j]]
    VectorNames = ['Disp_FEB']
    
    for i in range(len(VectorNames)) :
        arrayVector = vtk.util.numpy_support.numpy_to_vtk(VectorData[i], deep=True)
        arrayVector.SetName(VectorNames[i])
        dataVector = polydata.GetPointData()
        dataVector.AddArray(arrayVector)
        dataVector.Modified()
        
    fname = './NewFiles/test_' +str(j)+'.vtk'
    directory = os.path.dirname(fname)
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    writer = vtk.vtkDataSetWriter()
    writer.SetFileName(fname)
    writer.SetInputData(polydata)
    print('Writing ',fname)
    writer.Write()