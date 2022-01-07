
import numpy as np
import os
import sys
from anytree import NodeMixin, RenderTree
from anytree import LevelOrderIter
from inspect import currentframe, getframeinfo
from vtk.util.numpy_support import vtk_to_numpy
import vtk

frameinfo = getframeinfo(currentframe())

lookup =       {'00464542': 'FEB',
                '01000000': 'ROOT',
                '01010000': 'HEADER',
                '01010001': 'VERSION',
                '01010002': 'NODES',
                '01010003': 'MAX_FACET_NODES',
                '01010004': 'COMPRESSION',
                '01010005': 'AUTHOR',
                '01010006': 'SOFTWARE',
                '01020000': 'DICTIONARY',
                '01021000': 'GLOBAL_VAR',
                '01022000': 'MATERIAL_VAR',
                '01023000': 'NODESET_VAR',
                '01024000': 'DOMAIN_VAR',
                '01025000': 'SURFACE_VAR',
                '01020001': 'DICTIONARY_ITEM',
                '01020002': 'ITEM_TYPE',
                '01020003': 'ITEM_FORMAT',
                '01020004': 'ITEM_NAME',
                '01020005': 'ITEM_UNKNOWN',
                '01030000': 'MATERIALS',
                '01030001': 'MATERIAL',
                '01030002': 'MATERIAL_ID',
                '01030003': 'MATERIAL_NAME',
                '01040000': 'MESH',
                '01041000': 'NODE_SECTION',
                '01041100': 'NODE_HEADER',
                '01041101': 'NODES',
                '01041102': 'DIM',
                '01041103': 'NAME',
                '01041200': 'NODE_COORDS',
                '01042000': 'DOMAIN_SECTION',
                '01042100': 'DOMAIN',
                '01042101': 'DOMAIN_HEADER',
                '01042102': 'ELEM_TYPE',
                '01042103': 'PART_ID',
                '01032104': 'ELEMENTS',
                '01032105': 'NAME',
                '01042200': 'ELEMENT_LIST',
                '01042201': 'ELEMENT',
                '01043000': 'SURFACE_SECTION',
                '01043100': 'SURFACE',
                '01043101': 'SURFACE_HEADER',
                '01043102': 'SURFACE_ID',
                '01043103': 'FACES',
                '01043200': 'FACE_LIST',
                '01043201': 'FACE',
                '01044000': 'NODESET_SECTION',
                '01044100': 'NODESET',
                '01044101': 'NODESET_HEADER',
                '01044102': 'NODESET_ID',
                '01044103': 'NODESET_NAME',
                '01044104': 'NODESET_SIZE',
                '01044200': 'NODELIST',
                '01045000': 'PARTS_SECTION',
                '01045100': 'PART',
                '01045101': 'PART_ID',
                '01045102': 'PART_NAME',
                '02000000': 'STATE_SECTION',
                '02010000': 'STATE_HEADER',
                '02010002': 'STATE_TIME',
                '02010003': 'STATE_UNKNOWN',
                '02020000': 'STATE_DATA',
                '02020001': 'STATE_VAR',
                '02020002': 'VARIABLE_ID',
                '02020003': 'VARIABLE_DATA',
                '02020100': 'GLOBAL_DATA',
                '02020200': 'MATERIAL_DATA',
                '02020300': 'NODE_DATA',
                '02020400': 'DOMAIN_DATA',
                '02020500': 'SURFACE_DATA',
                '02030000': 'STATE_UNKNOWN2',
                }

def asHexString(x):
    return "{:08x}".format(x)

class xplt_node(NodeMixin):
    def __init__(self,name,size=None,read=False,parent=None,children=None):
        self.name = name
        self.size = size
        self.read = read
        self.parent = parent
        self.debug = False
        if children:
            self.children = children

    def set(self):
        self.read = True
        if self.parent is not None:
            self.parent.set()

    def read_file(self,fid):
        chunk = np.fromfile(fid,dtype=np.uint32,count=1)
        if not lookup[asHexString(chunk[0])]==self.name:
            if self.debug:
                print("ID not met",self.name,lookup[asHexString(chunk[0])])
                print("Assuming",self.name,"not present and continuing")
            fid.seek(-4,1)
            #return
            raise ValueError("ID not met",self.name,lookup[asHexString(chunk[0])])
        if self.debug:
            print('Reading',self.name)
        if self.name != 'FEB':
            self.size = np.fromfile(fid,dtype=np.uint32,count=1)[0]
        if self.debug:
            print('Reading',self.name,self.size)
        if self.read == True:
            if self.is_leaf:
                if self.size>0:
                    setattr(self, self.name, globals()["CALL_"+self.name](fid,self))
                else:
                    if self.debug:
                        print(self.name, 'has a size',self.size,'skipping reading it')
            else:
                for c in self.children:
                    try:
                        c.read_file(fid)
                    except:
                        c.parent=None #essentially delete the node
        else:
            fid.seek(self.size,1)
        

    def rank(self):
        i=0
        for node in LevelOrderIter(self.parent.parent,filter_=lambda n: n.depth==self.depth and n.name==self.name):
            if node is self:
                return i
            i += 1
        return i

def find_nstates(fid):
    nstate = 1
    feb.read_file(fid)
    while(True):
        try:
            state_section[nstate-1].read_file(fid)
            nstate += 1
            state_section.append(xplt_node('STATE_SECTION',parent=feb))
        except:
            # print("Number of states is",nstate)
            return nstate

def read1int32(fid):
    return np.fromfile(fid,dtype=np.uint32,count=1)[0]

def readchar(fid,n):
    strbits = np.fromfile(fid, dtype=np.int8, count=n)
    a_string = ''
    for i in strbits:
        if i>0:
            a_string += chr(i)
        else:
            return a_string
    return a_string

def CALL_DIM(fid,node):
    nDim = read1int32(fid)
    # print("Number of dimensions is",nDim)
    return nDim

def CALL_SOFTWARE(fid,node):
    chunk = np.fromfile(fid,dtype=np.int32,count=1)
    software = readchar(fid,chunk[0])
    # print('Software is',software)
    return software

def CALL_ITEM_TYPE(fid,node):
    itemtype = read1int32(fid)
    # print('Item type is',itemtype)
    return itemtype

def CALL_ELEMENTS(fid,node):
    nElem = read1int32(fid)
    # print('Number of elements is',nElem)
    return nElem

def CALL_ELEM_TYPE(fid,node):
    elementType = read1int32(fid)
    # print('Element type is',elementType)
    return elementType

def CALL_ITEM_FORMAT(fid,node):
    itemFormat = read1int32(fid)
    # print('Item format is',itemFormat)
    return itemFormat

def CALL_ITEM_UNKNOWN(fid,node):
    chunk = np.fromfile(fid,dtype=np.int32,count=1)
    itemUnknown = readchar(fid,chunk[0])
    # print('Unknown is',itemUnknown )
    return itemUnknown

def CALL_ITEM_NAME(fid,node):
    itemName = readchar(fid,64)
    # print('Item name is',itemName)
    return itemName

def CALL_NODES(fid,node):
    nNodes = read1int32(fid)
    # print("Number of nodes is",nNodes)
    return nNodes

def CALL_NODE_COORDS(fid,node):
    nNodes = node.path[0].children[1].children[0].children[0].children[0].NODES
    nDim = 3
    for i in range(nNodes):
        coords = np.fromfile(fid,dtype=np.float32,count=nDim)
    return coords

def CALL_STATE_TIME(fid,node):
    state_time = np.fromfile(fid,dtype=np.float32,count=1)[0]
    # print("State time is",state_time)
    return state_time

def CALL_VARIABLE_ID(fid,node):
    variableId = read1int32(fid)
    # print('Variable_ID is',variableId)
    return variableId

def CALL_VARIABLE_DATA(fid,node):
    VarId = node.parent.children[0].VARIABLE_ID-1
    
    if node.parent.parent.name == 'NODE_DATA':
        nPts = node.path[0].children[1].children[0].children[0].children[0].NODES
        ItemType = node.path[0].children[0].children[1].children[0].children[VarId].children[0].ITEM_TYPE
    elif node.parent.parent.name == 'DOMAIN_DATA':
        nPts = node.path[0].children[1].children[1].children[0].children[0].children[2].ELEMENTS
        ItemType = node.path[0].children[0].children[1].children[1].children[VarId].children[0].ITEM_TYPE
        
    if ItemType == 0:
        DataSize = 1
    elif ItemType == 1:
        DataSize = nPts*3+1
    elif ItemType == 2:
        DataSize = nPts*6+1
    elif ItemType == 3:
        DataSize = nPts*3
    elif ItemType == 4:
        DataSize = nPts*16
    elif ItemType == 5:
        DataSize = nPts*9
        
    # print(node.parent.parent.name, node.rank())
    chunk = np.fromfile(fid,dtype=np.uint32,count=1) #unsure what this refers to
    if node.parent.parent.name == 'NODE_DATA':
        chunk = np.fromfile(fid,dtype=np.float32,count=DataSize)
    elif node.parent.parent.name == 'DOMAIN_DATA':
        chunk = np.fromfile(fid,dtype=np.float32,count=DataSize)
    return chunk

    
def ReadObj(obj,Data,DataOut):
    # print(StateCount)
    if obj.name == 'FEB':
        if Data == 'VARIABLE_DATA':
            itemNames = ReadObj(obj,'ITEM_NAME',[])
            itemTypes = ReadObj(obj,'ITEM_TYPE',[])
            stateTimes = ReadObj(obj, 'STATE_TIME', [])
            DataOut.append(itemNames)
            DataOut.append(itemTypes)
            DataOut.append(stateTimes)
            for i in range(len(itemNames)):
                DataOut.append([])
                for j in range(len(stateTimes)):
                    DataOut[i+3].append([])
    if obj.read:
        if obj.is_leaf:
            if obj.name == Data:
                if Data == 'STATE_TIME':
                    # print('Time is ', eval('obj.'+Data))
                    DataOut.append(eval('obj.'+Data)) 
                elif Data == 'VARIABLE_DATA':
                    pathNames =[]
                    for i in range(len(obj.path)):
                        pathNames.append(obj.path[i].name)
                    # if 'NODE_DATA' in pathNames:
                    #     DataOut.append(eval('obj.'+Data)) 
                    nItems = len(DataOut[0])
                    nStates = len(DataOut[2])
                    for j in range(nStates):
                        for i in range(nItems):
                            if not DataOut[i+3][j]:
                                if hasattr(obj, Data):
                                    DataOut[i+3][j].append(eval('obj.'+Data)) 
                                else:
                                    DataOut[i+3][j].append('N/A')
                                break
                        else:
                            continue
                        break
                elif Data == 'ITEM_NAME':
                    DataOut.append(eval('obj.'+Data))
                elif Data == 'ITEM_TYPE':
                    DataOut.append(eval('obj.'+Data))
                elif Data == 'NODES':
                    DataOut.append(eval('obj.'+Data))
                elif Data == 'ELEMENTS':
                    DataOut.append(eval('obj.'+Data))
        else:
            for i in range(len(obj.children)):
                # print(StateCount)
                DataOut = ReadObj(obj.children[i],Data,DataOut)
           
    return DataOut
    
        
class xpltObj:
    def __init__(self,filename,filesize,nstates, feb):
        self.filename = filename
        self.filesize = filesize
        self.nstates  = nstates
        self.tree     = feb
        
        self.nNodes   = ReadObj(feb,'NODES',[])[0]
        self.nElems   = ReadObj(feb,'ELEMENTS',[])[0]
        self.VarData  = ReadObj(feb,'VARIABLE_DATA',[])
        self.VarNames = self.VarData[0]
        self.nVar     = len(self.VarData[0])
        self.VarType  = self.VarData[1]
        self.StateTimes = self.VarData[2]
        
        
        # Get State Times
        State_Times = np.zeros(nstates)
        for i in range(nstates):
            State_Times[i] = state_time[i].STATE_TIME
        
        
        # Get Variable 2
        self.Stress_X  = np.zeros((nstates,self.nElems))
        self.Stress_Y  = np.zeros((nstates,self.nElems))
        self.Stress_Z  = np.zeros((nstates,self.nElems))
        self.Stress_XY = np.zeros((nstates,self.nElems))
        self.Stress_YZ = np.zeros((nstates,self.nElems))
        self.Stress_XZ = np.zeros((nstates,self.nElems))
        self.Disp_X    = np.zeros((nstates,self.nNodes))
        self.Disp_Y    = np.zeros((nstates,self.nNodes))
        self.Disp_Z    = np.zeros((nstates,self.nNodes))
        
        for i in range(self.nstates):
            for j in range(self.nVar):
                if self.VarData[0][j] =='stress':
                    break
            stress = self.VarData[j+3][i][0]
            for j in range(self.nElems):
                self.Stress_X[i,j]  = stress[j*6+1]
                self.Stress_Y[i,j]  = stress[j*6+2]
                self.Stress_Z[i,j]  = stress[j*6+3]
                self.Stress_XY[i,j] = stress[j*6+4]
                self.Stress_YZ[i,j] = stress[j*6+5]
                self.Stress_XZ[i,j] = stress[j*6+6]
            for j in range(self.nVar):
                if self.VarData[0][j] =='displacement':
                    break
            displacement = self.VarData[j+3][i][0]
            for j in range(self.nElems):
                self.Disp_X[i,j]  = displacement[j*3+1]
                self.Disp_Y[i,j]  = displacement[j*3+2]
                self.Disp_Z[i,j]  = displacement[j*3+3]
        

if __name__ == '__main__':
    
    filename = sys.argv[1]
    
    feb = xplt_node('FEB')
    feb.set()
    root = xplt_node('ROOT',parent=feb)
    mesh = xplt_node('MESH',parent=feb)
    state_section = []
    state_section.append(xplt_node('STATE_SECTION',parent=feb))
    
    xpltFile = 'resample.xplt'
    fid = open(xpltFile, 'rb')
    fid.seek(0,2)
    file_size=fid.tell()
    fid.seek(0)
    nstates = find_nstates(fid)
    fid.seek(0)
    
    header = xplt_node('HEADER',parent=root)
    dictionary = xplt_node('DICTIONARY',parent=root)
    
    version = xplt_node('VERSION',parent=header)
    max_facet_nodes = xplt_node('MAX_FACET_NODES',parent=header)
    compression = xplt_node('COMPRESSION',parent=header)
    author = xplt_node('AUTHOR',parent=header)
    software = xplt_node('SOFTWARE',parent=header)
    
    global_var  = xplt_node('GLOBAL_VAR',parent=dictionary)
    nodeset_var  = xplt_node('NODESET_VAR',parent=dictionary)
    domain_var  = xplt_node('DOMAIN_VAR',parent=dictionary)
    surface_var  = xplt_node('SURFACE_VAR',parent=dictionary)
    
    dictionary_item_global = []
    dictionary_item_nodeset = []
    dictionary_item_domain = []
    dictionary_item_surface = []
    item_type = []
    item_format = []
    item_unknown = []
    item_name = []
    MAXVAR = 2
    for i in range(MAXVAR): #assume there are maximum MAXVAR dictionary items for each (global, nodeset, domain, and surface)
        dictionary_item_global.append(xplt_node('DICTIONARY_ITEM',parent=global_var))
        item_type.append(xplt_node('ITEM_TYPE',parent=dictionary_item_global[i]))
        item_format.append(xplt_node('ITEM_FORMAT',parent=dictionary_item_global[i]))
        item_unknown.append(xplt_node('ITEM_UNKNOWN',parent=dictionary_item_global[i]))
        item_name.append(xplt_node('ITEM_NAME',parent=dictionary_item_global[i]))
        item_type[-1].set()
        item_name[-1].set()
        item_format[-1].set()
    
        dictionary_item_nodeset.append(xplt_node('DICTIONARY_ITEM',parent=nodeset_var))
        item_type.append(xplt_node('ITEM_TYPE',parent=dictionary_item_nodeset[i]))
        item_format.append(xplt_node('ITEM_FORMAT',parent=dictionary_item_nodeset[i]))
        item_unknown.append(xplt_node('ITEM_UNKNOWN',parent=dictionary_item_nodeset[i]))
        item_name.append(xplt_node('ITEM_NAME',parent=dictionary_item_nodeset[i]))
        item_type[-1].set()
        item_name[-1].set()
        item_format[-1].set()
    
        dictionary_item_domain.append(xplt_node('DICTIONARY_ITEM',parent=domain_var))
        item_type.append(xplt_node('ITEM_TYPE',parent=dictionary_item_domain[i]))
        item_format.append(xplt_node('ITEM_FORMAT',parent=dictionary_item_domain[i]))
        item_unknown.append(xplt_node('ITEM_UNKNOWN',parent=dictionary_item_domain[i]))
        item_name.append(xplt_node('ITEM_NAME',parent=dictionary_item_domain[i]))
        item_type[-1].set()
        item_name[-1].set()
        item_format[-1].set()
    
        dictionary_item_surface.append(xplt_node('DICTIONARY_ITEM',parent=surface_var))
        item_type.append(xplt_node('ITEM_TYPE',parent=dictionary_item_surface[i]))
        item_format.append(xplt_node('ITEM_FORMAT',parent=dictionary_item_surface[i]))
        item_unknown.append(xplt_node('ITEM_UNKNOWN',parent=dictionary_item_surface[i]))
        item_name.append(xplt_node('ITEM_NAME',parent=dictionary_item_surface[i]))
        item_type[-1].set()
        item_format[-1].set()
    
    node_section    = xplt_node('NODE_SECTION',parent=mesh)
    domain_section  = xplt_node('DOMAIN_SECTION',parent=mesh)
    surface_section = xplt_node('SURFACE_SECTION',parent=mesh)
    nodeset_section = xplt_node('NODESET_SECTION',parent=mesh)
    parts_section   = xplt_node('PARTS_SECTION',parent=mesh)
    
    node_header = xplt_node('NODE_HEADER',parent=node_section)
    node_coords = xplt_node('NODE_COORDS',parent=node_section)
    
    domain = xplt_node('DOMAIN',parent=domain_section)
    
    nodes = xplt_node('NODES',parent=node_header)
    dim   = xplt_node('DIM',parent=node_header)
    name  = xplt_node('NAME',parent=node_header) 
    
    domain_header = xplt_node('DOMAIN_HEADER',parent=domain)
    element_list = xplt_node('ELEMENT_LIST',parent=domain)
    
    elem_type = xplt_node('ELEM_TYPE',parent=domain_header)
    part_id = xplt_node('PART_ID',parent=domain_header)
    elements = xplt_node('ELEMENTS',parent=domain_header)
    name = xplt_node('NAME',parent=domain_header) 
    
    state_header = []
    state_unknown2 = []
    state_data = []
    state_time = []
    state_unknown = []
    global_data = []
    node_data = []
    domain_data = []
    surface_data = []
    state_var_node = []
    variable_id_node = []
    variable_data_node = []
    state_var_domain = []
    variable_id_domain = []
    variable_data_domain = []
    
    for i in range(nstates):
        state_header.append(xplt_node('STATE_HEADER',parent=state_section[i]))
        state_unknown2.append(xplt_node('STATE_UNKNOWN2',parent=state_section[i]))
        state_data.append(xplt_node('STATE_DATA',parent=state_section[i]))
    
        state_time.append(xplt_node('STATE_TIME',parent=state_header[i]))
        state_unknown.append(xplt_node('STATE_UNKNOWN',parent=state_header[i]))
    
        global_data.append(xplt_node('GLOBAL_DATA',parent=state_data[i]))
        node_data.append(xplt_node('NODE_DATA',parent=state_data[i]))
        domain_data.append(xplt_node('DOMAIN_DATA',parent=state_data[i]))
        surface_data.append(xplt_node('SURFACE_DATA',parent=state_data[i]))
    
        for j in range(MAXVAR):
            state_var_node.append(xplt_node('STATE_VAR',parent=node_data[i]))
            variable_id_node.append(xplt_node('VARIABLE_ID',parent=state_var_node[-1]))
            variable_data_node.append(xplt_node('VARIABLE_DATA',parent=state_var_node[-1]))
            variable_data_node[-1].set()
            variable_id_node[-1].set()
    
            state_var_domain.append(xplt_node('STATE_VAR',parent=domain_data[i]))
            variable_id_domain.append(xplt_node('VARIABLE_ID',parent=state_var_domain[-1]))
            variable_data_domain.append(xplt_node('VARIABLE_DATA',parent=state_var_domain[-1]))
            variable_id_domain[-1].set()
        variable_data_domain[-2].set()
        variable_data_domain[-1].set()
    
    nodeset_var.set()
    domain_var.set()
    software.set()
    nodes.set()
    dim.set()
    domain.set()
    # node_coords.set()
    elem_type.set()
    elements.set()
    
    for i in range(nstates):
        state_time[i].set()
    
    feb.read_file(fid)
    
    # for pre,fill,node in RenderTree(feb):
    #     treestr = u"%s%s" %(pre,node.name)
    #     print(treestr.ljust(8),node.name,node.read)
    
    
    xplt_obj = xpltObj(filename,file_size,nstates, feb)
    
    print('Filename: ', xplt_obj.filename)
    print('File size: ', xplt_obj.filesize)
    print('Number of States: ', xplt_obj.nstates)
    print('The tree of data is called: ', 'feb')
    print('Number of Nodes: ', xplt_obj.nNodes)
    print('Number of Elements: ', xplt_obj.nElems)
    print('The saved variables are ', xplt_obj.VarNames)
        
    reader = vtk.vtkPolyDataReader()
    
    # Read the source file.
    reader.SetFileName('resample.vtk')
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    
    polydata = reader.GetOutput()
    
    # Add Vector Data on Points
    VectorData = [xplt_obj.Stress_X[0,:],
                  xplt_obj.Stress_Y[0,:],
                  xplt_obj.Stress_Z[0,:],
                  xplt_obj.Stress_XY[0,:],
                  xplt_obj.Stress_YZ[0,:],
                  xplt_obj.Stress_XZ[0,:]]
    
    VectorNames = ['Stress_X','Stress_Y','Stress_Z','Stress_XY','Stress_YZ','Stress_XZ']
    
    for i in range(len(VectorNames)) :
        arrayVector = vtk.util.numpy_support.numpy_to_vtk(VectorData[i], deep=True)
        arrayVector.SetName(VectorNames[i])
        dataVectors = polydata.GetCellData()
        dataVectors.AddArray(arrayVector)
        dataVectors.Modified()
        
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
    PointData = [Stress_X_Pt,Stress_Y_Pt,Stress_Z_Pt,Stress_XY_Pt,Stress_YZ_Pt,Stress_XZ_Pt]
    PointNames = ['Stress_X_Pt','Stress_Y_Pt','Stress_Z_Pt','Stress_XY_Pt','Stress_YZ_Pt','Stress_XZ_Pt']
    for i in range(len(PointNames)) :
        arrayPoint = vtk.util.numpy_support.numpy_to_vtk(PointData[i], deep=True)
        arrayPoint.SetName(PointNames[i])
        dataPoints = polydata.GetPointData()
        dataPoints.AddArray(arrayPoint)
        dataPoints.Modified()
                
    fname = './test.vtk'
    directory = os.path.dirname(fname)
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    writer = vtk.vtkDataSetWriter()
    
    writer.SetFileName(fname)
    writer.SetInputData(polydata)
    print('Writing ',fname)
    writer.Write()
    
    
    