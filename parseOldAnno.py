import xml.etree.ElementTree as ElTre
import os
from sys import argv
import json


def main(inpath, outpath):
    proteome = {}
    count = {}
    clan_dict = {}
    tools = ['pfam', 'smart', 'tmhmm', 'coils2', 'flps', 'signalp', 'seg']
    for tool in tools:
        if tool == 'coils2':
            proteome, count, clan_dict = xmlreader(inpath + '/coils.xml', tool, proteome, clan_dict, count)
        else:
            proteome, count, clan_dict = xmlreader(inpath + '/' + tool + '.xml', tool, proteome, clan_dict, count)
    for protein in proteome:
        for tool in tools:
            if tool not in proteome[protein]:
                proteome[protein][tool] = {}
    outdict = {'feature': proteome, 'count': count, 'clan': clan_dict}
    save2json(outdict, outpath)


def save2json(dict2save, outpath):
    jsonOut = json.dumps(dict2save, ensure_ascii=False)
    f = open(outpath, 'w')
    f.write(jsonOut)
    f.close()


def xmlreader(path, tool, proteome, clan_dict, count):
    start = 0
    if os.path.exists(path):
        try:
            xmltree = ElTre.parse(path)
        except ElTre.ParseError:
            raise Exception(path + ' is not a legit xml file')
        root = xmltree.getroot()
        for protein in root:
            p_id = protein.attrib["id"]
            plength = protein.attrib["length"]
            # set up of protein IDs, differentiate between proteins from different files
            if not (p_id in proteome):
                proteome[p_id] = {'length': int(plength), tool: {}}
            else:
                proteome[p_id][tool] = {}
            # set up of datastructure to store annotations

            for feature in protein:
                if len(feature) > 0:
                    ftype = tool + "_" + feature.attrib["type"].replace(' ', '_')
                    feat_eval = 'NA'
                    if ftype == 'signalp_SIGNAL':
                        ftype = 'signalp_SIGNALP'
                    elif ftype == 'coils2_coiled_coil':
                        ftype = 'coils_coiled_coil'

                    # evalue check: family/ftype based
                    if ('evalue' in feature.attrib and float(feature.attrib["evalue"]) <= 0.001) or 'evalue' not in feature.attrib:
                        if 'evalue' in feature.attrib:
                            feat_eval = float(feature.attrib["evalue"])
                        if 'clan' in feature.attrib and ftype not in clan_dict and not feature.attrib["clan"] == '---':
                            fclan = feature.attrib["clan"]
                            clan_dict[ftype] = fclan
                        proteome[p_id][tool][ftype] = {'instance': [], 'evalue': feat_eval}

                        i = 0
                        # counting appended instances
                        inst_count = 0
                        for instance in feature:
                            inst_eval = 'NULL'
                            # XMLcase 1: feature instance contains evalue information (XML field inst_eval)
                            if 'inst_eval' in instance.attrib:
                                # print tool + " instance evalue: "+ str(instance.attrib)
                                inst_eval = float(instance.attrib["inst_eval"])
                                start = int(instance.attrib["start"])
                                end = int(instance.attrib["end"])

                                if inst_eval <= 0.01:
                                    proteome[p_id][tool][ftype]['instance'].append((start, end, inst_eval))
                                    inst_count += 1

                            # XMLcase 2: feature instance contains NO evalue information (XML field inst_eval)
                            else:
                                # NO instance based evalue information --> no instances can be rejected:
                                # set inst_count = 1
                                if len(instance.attrib) == 2:
                                    start = int(instance.attrib["start"])
                                    end = int(instance.attrib["end"])
                                    proteome[p_id][tool][ftype]['instance'].append((start, end, 'NA'))
                                    inst_count += 1
                                else:
                                    if i == 0:
                                        start = int(instance.attrib["start"])
                                        i = 1
                                        inst_count += 1
                                    else:
                                        end = int(instance.attrib["end"])
                                        proteome[p_id][tool][ftype]['instance'].append((start, end, 'NA'))
                                        i = 0
                        # any instance appended?
                        if inst_count < 1:
                            proteome[p_id][tool].pop(ftype)
                        elif ftype in count:
                            count[ftype] += inst_count
                        else:
                            count[ftype] = inst_count
    else:
        raise Exception(path + " does not exist")
    return proteome, count, clan_dict


if __name__ == '__main__':
    main(argv[1], argv[2])
