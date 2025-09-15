import sys
import vtk
import glob
import re


# member variables are dynamically added by parsing data files
# for Sylinder, Protein, and ConBlock classes

class Sylinder(object):
    end0 = None
    end1 = None
    pass


class ConBlock(object):
    end0 = None
    end1 = None
    pass


class Frame:

    def __init__(self, sylinderFile=None,
                 conBlockFile=None):
        self.sylinders = []
        self.conBlocks = []
        self.parseSylinderFile(sylinderFile)
        self.parseConBlockFile(conBlockFile)

    def parseFile(self, dataFile, objType, objList):
        print("Parsing data from " + dataFile)
        # create vtk reader
        reader = vtk.vtkXMLPPolyDataReader()
        reader.SetFileName(dataFile)
        reader.Update()
        data = reader.GetOutput()

        # fill data
        # step 1, end coordinates
        nObj = int(data.GetPoints().GetNumberOfPoints() / 2)
        print("parsing data for ", nObj, " sylinders")
        for i in range(nObj):
            s = objType()
            s.end0 = data.GetPoints().GetPoint(2 * i)
            s.end1 = data.GetPoints().GetPoint(2 * i + 1)
            objList.append(s)

        # step 2, member cell data
        numCellData = data.GetCellData().GetNumberOfArrays()
        print("Number of CellDataArrays: ", numCellData)
        for i in range(numCellData):
            cdata = data.GetCellData().GetArray(i)
            dataName = cdata.GetName()
            print("Parsing Cell Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName, cdata.GetTuple(j))

        # step 3, member point data
        numPointData = data.GetPointData().GetNumberOfArrays()
        print("Number of PointDataArrays: ", numPointData)
        for i in range(numPointData):
            pdata = data.GetPointData().GetArray(i)
            dataName = pdata.GetName()
            print("Parsing Point Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName + "0", pdata.GetTuple(2 * j))
                setattr(objList[j], dataName + "1", pdata.GetTuple(2 * j + 1))

        # output all data for debug
        for s in objList[:10]:
            # print(s.end0, s.end1)
            attrs = vars(s)
            print('*************************************')
            print('\n'.join("%s: %s" % item for item in attrs.items()))
            print('*************************************')

        print("-------------------------------------")

    def parseSylinderFile(self, sylinderFile):
        self.parseFile(sylinderFile, Sylinder, self.sylinders)

    def parseConBlockFile(self, conBlockFile):
        self.parseFile(conBlockFile, ConBlock, self.conBlocks)


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


def doSomethingPerFrame(frame):
    pass


def main():
    # get file list
    SylinderFileList = glob.glob('./result/result*/Sylinder_*.pvtp')
    ConBlockFileList = glob.glob('./result/result*/ConBlock_*.pvtp')

    # sort as numerical order
    SylinderFileList.sort(key=getFrameNumber_lambda)
    ConBlockFileList.sort(key=getFrameNumber_lambda)

    print(SylinderFileList)
    print(ConBlockFileList)

    assert len(SylinderFileList) == len(ConBlockFileList)

    # example
    for i in range(len(SylinderFileList)):
        # get frame
        frame = Frame(SylinderFileList[i], ConBlockFileList[i])
        # do something
        doSomethingPerFrame(frame)


if __name__ == '__main__':
    main()
