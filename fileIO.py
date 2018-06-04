import copy
import numpy as np
#from .. element import Element


class LatticeFile(object):
    '''
    It can be reused by all supported external files
    elementList: List of all defined elements
    elementNameDict: Look up table of elements from its name {'name': index}
    beamlineList: List of all defined beamlines
    beamlineNameDict: Look up table of beamlines from its name

    '''
    def __init__(self):
        self.lattice_format='PyLatte'
        self.elementList = []
        self.elementNameDict = {}
        self.beamlineList = []
        self.beamlineNameDict = {}
        self.beamlineElementSet = {}
        self.useLineList = []
        self.elementPosInUseLine = None

        self.useline= ''

    def checkType(self, typename, parameterName=None):
        '''
        Check if the type is valid, should be inherit by all child for there own definition.
        This parent will adopt the definition in .. element
        :param typename:
        :param parameterName:
        :return: True if type and paramterName is allowed, otherwise False
        '''
        if typename.upper() not in Element.elementTypes:
            return False
        else:
            if parameterName is None:
                return True
            else:
                return parameterName.upper() in Element.get_property_name(typename.upper())

    def toLatte(self):
        pass
    def fromLatte(self):
        pass
    def parseFrom(self, filename):
        pass

    def getElementIndex(self, elename):
        '''
        From the name to get the index in the list of elementNameDict;  Ignore the leading '-'

        :param elename:
        :return: the index in elementNameDict or None if name not found
        '''
        elename = elename.upper()
        if elename[0] == '-':
            return self.getElementIndex(elename[1:])
        elif elename in self.elementNameDict:
            return self.elementNameDict[elename]
        else:
            return None


    def addElement(self, _name, _eletype, **params):
        '''
        Add a new element.
        :param _name: Name of the element to be added
        :param _eletype: Type of the element to be added, either checked by checkType(_eletype) or existing element name.
        :param params: A dict input on the parameterlist, The allowed parameter name checked by checkType(_eletype, param_name)
        :return: None
        '''
        roottype = ''
        #if self.checkType(_eletype) == False:
        #    if (self.getElementIndex(_eletype) is None):
        #        print ("The element {} with type {} is not recognized when adding this element".format(_name, _eletype))
        #        raise KeyError
        #    else:
        #        roottype = self.getElementRootType(_eletype)
        if self.getElementIndex(_eletype) is None:
            if self.checkType(_eletype) == False:
                print ("The element {} with type {} is not recognized when adding this element".format(_name, _eletype))
                raise KeyError
        else:
            roottype = self.getElementRootType(_eletype)


        thisele = {}
        name = _name.upper()
        eletype = _eletype.upper()
        thisele['NAME'] = name
        thisele['TYPE'] = eletype
        thisele['__USAGE_DICT'] = {}
        thisele['__ID_IN_USELINE'] = []

        for k, v in params.items():
            if k.upper() == 'NAME':
                pass
            elif k.upper() == 'TYPE':
                pass
            elif self.checkType(eletype, k) or self.checkType(roottype, k):
                thisele[k.upper()] = v
            else:
                print('Unrecognized parameter name {} in element {} when adding'.format(k, name))
                raise KeyError

        ind = self.getElementIndex(name)
        if ind is None:
            cur_len = len(self.elementList)
            self.elementNameDict[name] = cur_len
            self.elementList.append(thisele)
        else:
            self.elementList[ind] = thisele
            print ("Warning, the element {} is redefined when adding this element, overwriting with newer one.".format(name))
        return

    def getElementRootType(self, elename):
        '''
        Get the type of the element even if it is derived from other files.
        :param elename: Name of the element
        :return: The root type of the element, if the type of the element is another element.
        '''
        elename = elename.upper()
        ind = self.getElementIndex(elename)
        if ind is None:
            print('Unrecognized element name {} when getting root type'.format(elename))
            raise KeyError
        cur_type = self.elementList[ind]['TYPE']
        if cur_type in self.elementNameDict:
            return self.getElementRootType(cur_type)
        else:
            return cur_type


    def getElementProperties(self, elename, keyname=None, partial_key=''):
        '''
        Return the element's property value based on its parameter name or partial name.
        :param elename: The name of the element
        :param keyname: The parameter name, if empty, all properties,filtered by the partial_key will be returned
        :param partial_key: Only the parameter name contains partial_key will be returned.
        :return: The parameter value, if keyname is not None, or dict of {parameterName: parameterValue} is returned
        '''
        elename = elename.upper()
        ind = self.getElementIndex(elename)
        if ind is None:
            print('Unrecognized element name {} when getting properties'.format(elename))
            raise KeyError
        if keyname is None:
            root_type = self.getElementRootType(elename)
            cur_type = self.elementList[ind]['TYPE']

            if cur_type == root_type:
                return {k: v for k, v in self.elementList[ind].items()
                        if partial_key.upper() in k and k not in ['NAME', '__USAGE_DICT', '__ID_IN_USELINE']}
            else:
                temp = self.getElementProperties(cur_type, None, partial_key)
                temp.update({k: v for k, v in self.elementList[ind].items()
                     if partial_key.upper() in k and k not in ['TYPE', 'NAME', '__USAGE_DICT', '__ID_IN_USELINE']})
                return temp

        else:
            keyname = keyname.upper()
            if keyname in self.elementList[ind]:
                return self.elementList[ind][keyname]
            elif self.elementList[ind]['TYPE'] in self.elementNameDict:
                return self.getElementProperties(self.elementList[ind]['TYPE'], keyname)
            else:
                return None



    def getParentElements(self, elename):
        pe_list = []
        ind = self.getElementIndex(elename)
        if ind is None:
            print('Unrecognized element name {} when finding parents'.format(elename))
            raise KeyError
        t = self.elementList[ind]['TYPE']
        pe_list.append(t)
        if t in self.elementNameDict:
            pe_list += self.getParentElements(t)
        return pe_list

    def compareElements(self, other_lattice, elename, other_name=None):
        elename = elename.upper()
        if other_name is None:
            other_name = elename
        else:
            other_name = other_name.upper()
        ind_other = other_lattice.getElementIndex(other_name)
        ind = self.getElementIndex(elename)
        if ind is None or ind_other is None:
            return False
        return self.getElementProperties(elename) == other_lattice.getElementProperties(other_name)

    def modifyElement(self, elename, increment=False, **params):
        '''
        Modify properties of one element
        :param elename: The element name to be modified
        :param increment: Flag to choose overwrite (False, default) or add to existing value (True)
        :param params: dictionary contains the changes.
        :return: None
        '''
        elename = elename.upper()
        ind = self.getElementIndex(elename)
        if ind is None:
            print('Unrecognized element name {} when modifying'.format(elename))
            raise KeyError
        for k, v in params.items():
            if self.checkType(self.elementList[ind]['TYPE'], k):
                if k in self.elementList[ind] and increment:
                    self.elementList[ind][k.upper()] += v
                else:
                    self.elementList[ind][k.upper()] = v
            else:
                print('Unrecognized parameter name {} in element {}'.format(k, elename))
                raise KeyError
        return

    def modifyAllElements(self, eletype, increment=False, name_contain='', **params):
        eletype = eletype.upper()

        for ele in self.elementList:
            if ele['TYPE']==eletype and name_contain.upper() in ele['NAME']:
                self.modifyElement(ele['NAME'], increment, **params)

                # parameterName = parameterName.upper()
                # if eletype in self.elementNameDict:
                #    roottype = self.getElementRootType(eletype)
                # else:
                #   roottype = eletype
                # if self.checkType(eletype):
                #    roottype = eletype
                # else:
                #    roottype = self.getElementRootType(eletype)
                # if not self.checkType(roottype, parameterName):
                #    print('The parameter name {} is not good for element type {}'.format(parameterName, eletype))
                #    raise KeyError
                #current_value = self.getElementProperties(ele['NAME'], parameterName)
                #future_value = parameterValue
                #if increment and current_value is not None:
                #    future_value += current_value
                #if ele['TYPE'] == eletype:
                #    ele[parameterName] = future_value
                #elif eletype in self.getParentElements(ele['NAME']) and current_value != future_value:
                #    ele[parameterName] = future_value

        #for ele in self.elementList:
        #    if eletype in self.getParentElements(ele['NAME']) and name_contain in ele['NAME'] and ele['TYPE']!=eletype \
        #            and self.getElementProperties(ele['NAME'],parameterName) != parameterValue:
        #        ele[parameterName] = parameterValue
            #elif self.elementList[i]['TYPE'] in self.elementNameDict:
            #    if self.elementList[self.elementNameDict[self.elementList[i]['TYPE']]] == eletype:
            #        self.elementList[i][parameterName]=parameterValue

    def isDrift(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'DRIFT' in parent_type:
            temp=self.getElementProperties(ele_name)
            return temp
        else:
            return False

    def isDipole(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'DIPOLE' in parent_type:
            temp=self.getElementProperties(ele_name)
            if 'K1' not in temp:
                temp['K1']=0
            return temp
        else:
            return False

    def isQuadrupole(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'QUADRUPOLE' in parent_type:
            temp=self.getElementProperties(ele_name)
            if 'K1' not in temp:
                temp['K1']=0
            return temp
        else:
            return False
    def isSolenoid(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'SOLENOID' in parent_type:
            temp=self.getElementProperties(ele_name)
            return temp
        else:
            return False
    def isCavity(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'CAVITY' in parent_type:
            temp=self.getElementProperties(ele_name)
            return temp
        else:
            return False

    def plotBeamline(self, plt_axis, colors=('DarkOrchid', 'Maroon', 'DeepSkyBlue', 'ForestGreen'),
                     heights=(1.2, 0.8, 0.5, 0.5), s_start=0, other_components = {}):
        '''
        Plot the beamline using shapes.
        :param plt_axis: matplotlib axis variable
        :param beamline_name: name of beamline to be plotted.
        :param colors: The color for shapes that represent cavity, dipole, Quad and solenoid
        :param heights: The height for shapes that represent cavity, dipole, Quad and solenoid
        :param s_start: offset of the starting point, for display only
        :param other_components: Dictionary for other components to plot.
        Format:{components:{'real_length':True, 'type':'on_axis', 'height':0.2, 'color':'b'}}
        :return:
        '''
        if self.elementPosInUseLine is None:
            self.setUseLine()
        bl_pos, bl_list = self.elementPosInUseLine, self.useLineList
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        from matplotlib.patches import Ellipse

        plt_axis.set_xlim([s_start, s_start + bl_pos[-1]])

        plt_axis.set_ylim([-1, 1])
        plt_axis.set_yticks([])
        # plt_axis.xaxis.set_ticks_position('none')
        plt_axis.axhline(0, color='black', lw=0.5)

        if other_components is not None:
            other_components = {k.upper(): v for k, v in other_components.items()}
        for i in range(len(bl_list)):
            start = bl_pos[i] + s_start
            ele_name = bl_list[i]
            #ele=self.elementList[self.elementNameDict[ele_name]]
            #typelist = self.getParentElements(ele_name)

            l_ele = bl_pos[i + 1] - bl_pos[i]
            #lasttype = self.elementList[indlist[-1]]['TYPE']
            #lasttype=typelist[-1]
            if l_ele > 0:
                tempdict = self.isCavity(ele_name)
                if tempdict:
                    #plt_axis.add_patch(
                    #    Rectangle((start, -heights[0] / 2.0), tempdict['L'], heights[0], ec=colors[0], fc=colors[0]))
                    plt_axis.add_patch(
                        Ellipse((start+tempdict['L']/2.0, 0), tempdict['L'], heights[0] ,linewidth=1, color=colors[0]))
                    continue
                tempdict = self.isDipole(ele_name)
                if tempdict:
                    shift = heights[1] * 0.2 * np.sign(tempdict['K1'])
                    plt_axis.add_patch(
                        Rectangle((start, -heights[1] / 2.0 + shift), tempdict['L'], heights[1], angle=0.0, ec=colors[1],
                                  fc='none'))
                    continue
                tempdict = self.isQuadrupole(ele_name)
                if tempdict:
                    shift = heights[2] * 0.5 * np.sign(tempdict['K1'])
                    plt_axis.add_patch(
                        Rectangle((start, -heights[2] / 2.0 + shift), tempdict['L'], heights[2], angle=0.0, ec=colors[2],
                                  fc='none'))
                    continue
                tempdict = self.isSolenoid(ele_name)
                if tempdict:
                    plt_axis.add_patch(
                        Rectangle((start, -heights[3] / 2.0), tempdict['L'], heights[3], angle=0.0, ec=colors[3], fc='none'))
                    continue

            elif not self.isDrift(ele_name):
                lasttype = self.getParentElements(ele_name)[-1]
                if lasttype in other_components:
                    p = other_components[lasttype]

                    if p.get('real_length', False):
                        pass
                    else:
                        if p.get('type', 'on_axis') == 'on_axis':
                            h = p.get('height', 0.4)
                            c = p.get('color', 'y')
                            plt_axis.axvline(start + l_ele / 2.0, ymin=0.5 - h / 4.0, ymax=0.5 + h / 4.0, color=c)
                        else:
                            h = p.get('height', 0.3)
                            c = p.get('color', 'k')
                            plt_axis.axvline(start + l_ele / 2.0, ymin=1.0-h/2, ymax=1.0, color=c)
                            plt_axis.axvline(start + l_ele / 2.0, ymin=0 / 2, ymax=h/2, color=c)





    def getBeamlineIndex(self, linename):
        linename = linename.upper()
        if linename[0] == '-':
            return self.getBeamlineIndex(linename[1:])
        elif linename in self.beamlineNameDict:
            return self.beamlineNameDict[linename]
        else:
            return None

    def appendToBeamline(self, linename, *elenames):
        '''
        Append elements to a line
        :param linename: The line name to be appended. If the line does not exist, it will be created first
        :param elenames: The list of elements to be used.
        :return: None
        '''
        linename = linename.upper()

        ind = self.getBeamlineIndex(linename)
        if ind is None:
            addaline = {}
            addaline['NAME'] = linename
            addaline['LINE'] = []
            addaline['__USAGE_DICT']={}
            cur_len = len(self.beamlineList)
            self.beamlineList.append(addaline)
            self.beamlineNameDict[linename] = cur_len
            ind=cur_len

        for elename in elenames:
            ind_ele = self.getElementIndex(elename)  # element is really an element
            ind_line = self.getBeamlineIndex(elename)  # element is a line?
            if ind_ele is None and ind_line is None:
                print('No element or beamline named {} are defined yet in line {} when appending to a line.'.format(elename, linename))
                raise KeyError
            else:
                self.beamlineList[ind]['LINE'].append(elename.upper())
                cur_pos=len(self.beamlineList[ind]['LINE'])-1
                if ind_line is not None:
                    # self.beamlineNameDict[linename]=checkline
                    # self.beamlineNameDict[elenames]=ind
                    # self.beamlineList[ind],self.beamlineList[checkline]=self.beamlineList[checkline],self.beamlineList[ind]
                    if ind_line > ind:
                        temp = self.beamlineList.pop(ind)
                        self.beamlineList.append(temp)
                        for k, v in self.beamlineNameDict.items():
                            if v > ind:
                                self.beamlineNameDict[k] = v - 1
                        self.beamlineNameDict[linename] = len(self.beamlineList) - 1
                        ind_line = self.getBeamlineIndex(elename)
                        ind = len(self.beamlineList) - 1

                    if linename not in self.beamlineList[ind_line]['__USAGE_DICT']:
                        self.beamlineList[ind_line]['__USAGE_DICT'][linename]=[]
                    self.beamlineList[ind_line]['__USAGE_DICT'][linename].append(cur_pos)
                if ind_ele is not None:
                    if linename not in self.elementList[ind_ele]['__USAGE_DICT']:
                        self.elementList[ind_ele]['__USAGE_DICT'][linename]=[]
                    self.elementList[ind_ele]['__USAGE_DICT'][linename].append(cur_pos)



    def elementInLine(self, elename, linename):
        '''
        check if the element with name elename is in the line of linename
        '''
        ind_ele=self.getElementIndex(elename)
        ind_line=self.getBeamlineIndex(elename)
        if ind_ele is not None:
            ind = ind_ele
            dict_temp = self.elementList[ind]['__USAGE_DICT']
        elif ind_line is not None:
            ind = ind_line
            dict_temp = self.beamlineList[ind]['__USAGE_DICT']
        else:
            return False
        if linename.upper() in dict_temp:
            return True
        else:
            for k,v in dict_temp.items():
                if self.getBeamlineIndex(k) is not None:
                    return self.elementInLine(k, linename)

    def addReverseLine(self, newlinename, linename):
        '''
        Add reversed line with newlinename
        '''
        count = self.getBeamlineIndex(linename)
        if count is None:
            print('The line with name {} can not be found, no reverse line can be added'.format(linename))
        self.beamlineNameDict[newlinename] = len(self.beamlineList)
        self.beamlineList.append(copy.deepcopy(self.beamlineList[count]))
        self.beamlineList[-1]['NAME'] = newlinename
        self.beamlineList[-1]['LINE'].reverse()
        for i in range(len(self.beamlineList[-1]['LINE'])):
            elename = self.beamlineList[-1]['LINE'][i]
            if elename[0] == '-':
                self.beamlineList[-1]['LINE'][i] = elename[1:]
            #elif self.checkAsymmetryType(self.elementList[self.elementNameDict[elename]]):
            #    self.beamlineList[-1]['LINE'][i] = '-' + elename
            #else:
            #    pass

    def _expandLine(self, line_ind):

        expandedLine=[]
        briefline = self.beamlineList[line_ind]['LINE']
        for elename in briefline:
            if elename in self.beamlineNameDict:
                expandedLine += self._expandLine(self.beamlineNameDict[elename])
            else:
                expandedLine.append(elename)
        return expandedLine

    def setUseLine(self, linename=None):
        if linename is None:
            linename = self.beamlineList[-1]['NAME']
            line_ind = len(self.beamlineList)-1
        else:
            linename=linename.upper()
            line_ind = self.getBeamlineIndex(linename)
        self.useline=linename

        self.elementPosInUseLine = [0.0, ]

        if line_ind is None:
            print('The beamline {} does no exist, can not prepare the line to be used'.format(linename))
            raise KeyError

        self.useLineList = self._expandLine(line_ind)

        ele_ind=0
        for ele in self.useLineList:
            ind_ele = self.getElementIndex(ele)
            self.elementList[ind_ele]['__ID_IN_USELINE'].append(ele_ind)
            l_ele = self.getElementProperties(ele, 'L')
            if l_ele==None:
                l_ele=0.0
            self.elementPosInUseLine.append(self.elementPosInUseLine[-1] + l_ele)
            ele_ind += 1


    def loadAnElement(self, fromlattice, elename, prefix='', suffix=''):
        elename = elename.upper()
        prefix = prefix.upper()
        suffix = suffix.upper()
        ind = fromlattice.getElementIndex(elename)
        if ind is not None:
            ele = copy.deepcopy(fromlattice.elementList[ind])
            ind_this = self.getElementIndex(prefix + elename + suffix)
            if ind_this is not None:
                compare=self.compareElements(fromlattice, prefix + ele['NAME']+ suffix, other_name=ele['NAME'])
                if compare:
                    return
                else:
                    print('Warning, the element {} has different definition'.format(prefix+elename+suffix))
                    return

            if self.checkType(ele['TYPE']):
                ele['NAME'] = prefix + ele['NAME']+ suffix
                ele['__USAGE_DICT'] = {}
                ele['__ID_IN_USELINE'] = []
                self.elementNameDict[ele['NAME']] = len(self.elementList)
                self.elementList.append(ele)
            elif self.getElementIndex(prefix + ele['TYPE']) is None:
                self.loadAnElement(fromlattice, ele['TYPE'], prefix, suffix)
                ele['NAME'] = prefix + ele['NAME']+ suffix
                ele['TYPE'] = prefix + ele['TYPE']+ suffix
                ele['__USAGE_DICT'] = {}
                ele['__ID_IN_USELINE'] = []
                self.elementNameDict[ele['NAME']] = len(self.elementList)
                self.elementList.append(ele)
            elif self.compareElements(fromlattice, prefix + ele['TYPE']+ suffix, other_name=ele['TYPE']):
                ele['NAME'] = prefix + ele['NAME'] + suffix
                ele['TYPE'] = prefix + ele['TYPE'] + suffix
                ele['__USAGE_DICT'] = {}
                ele['__ID_IN_USELINE'] = []
                self.elementNameDict[ele['NAME']] = len(self.elementList)
                self.elementList.append(ele)
                pass
            else:
                print('Warning, the root element {} has different definition'.format(ele['TYPE']))
                exit(-1)
                '''self.loadAnElement(fromlattice, ele['TYPE'], prefix)
                ele['NAME'] = prefix + ele['NAME']
                ele['TYPE'] = prefix + ele['TYPE']
                self.elementNameDict[ele['NAME']] = len(self.elementList)
                self.elementList.append(ele)'''


        else:
            print('Can not load element {} from lattice'.format(elename))
            exit(-1)

    def loadALine(self, fromlattice, linename, reverse=False, prefix='', newname='', suffix=''):
        linename = linename.upper()
        if newname == '':
            newname = linename
        else:
            newname = newname.upper()
        prefix = prefix.upper()
        suffix = suffix.upper()
        ind = fromlattice.getBeamlineIndex(linename)
        ind_this = self.getBeamlineIndex(prefix + newname+ suffix)
        templine = []
        if ind is not None and ind_this is None:
            theline = copy.deepcopy(fromlattice.beamlineList[ind]['LINE'])

            if reverse:

                theline.reverse()
                for i in range(len(theline)):
                    if theline[i][0] == '-':
                        theline[i] = theline[i][1:]
                    else:
                        pass
                        # print(theline)
                        # for elename in reversed(fromlattice.beamlineList[pos]['LINE']):
                        #    if fromlattice.checkElementName(elename) is not None:
                        #        self.loadAnElement(fromlattice, elename, prefix)
                        #        #self.addToBeamLine(prefix+newname, prefix+elename)
                        #        templine.append(prefix+elename)
                        #    elif fromlattice.checkBeamlineName(elename)[0]>0:
                        #        self.loadALine(fromlattice, elename, reverse=True, prefix=prefix)
                        # self.addToBeamLine(prefix+newname, prefix+elename)
                        #        templine.append(prefix+elename)

            for elename in theline:
                if fromlattice.getElementIndex(elename) is not None:
                    if elename[0] == '-':
                        self.loadAnElement(fromlattice, elename[1:], prefix, suffix)
                        templine.append(prefix + elename[1:]+ suffix)

                    else:

                        self.loadAnElement(fromlattice, elename, prefix, suffix)
                        templine.append(prefix + elename+ suffix)


                elif fromlattice.getBeamlineIndex(elename) is not None:

                    if elename[0] == '-':
                        self.loadALine(fromlattice, elename[1:], reverse=False, prefix=prefix, suffix=suffix)
                        templine.append('-' + prefix + elename[1:]+ suffix)
                    else:
                        self.loadALine(fromlattice, elename, reverse=False, prefix=prefix, suffix= suffix)
                        templine.append(prefix + elename+ suffix)
            self.appendToBeamline(prefix + newname+ suffix, *templine)
        elif ind is None:
                print("The line {} doesnot exist in the source".format(linename))
                exit()
        elif ind_this is not None:
                print("The line {} is already in the target".format(newname))
        else:
                print("Something weird happend in loadALine when loading {}".format(linename))
                exit()
        return prefix + newname + suffix

