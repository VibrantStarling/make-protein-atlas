#Convert ".ptm.pep.xml" TPP PTMprophet results output to .csv file for downstream analysis
#output="interact.ptm.pep.csv"


import xml.sax
import sys
import os

sep = "\t"
PEPTIDE_PROPHET_THRESHOLD = 0

#TODO
#Have a running log, count how many spectra were reported in the file, append to a dataset log file


#convert xml output file from PTM prophet to csv
def convert(input,file_prefix):
    class TPPHandler(xml.sax.ContentHandler):

        #global psm_dict


        def __init__(self):
            self.CurrentData = ""
            self.psm_dict = {}


        # Call when an element starts
        def startElement(self, tag, attributes):

            if tag == "spectrum_query":
                self.psm_dict = {}
                self.psm_dict["spectrum"] = attributes["spectrum"]
                self.psm_dict["rt"] = attributes["retention_time_sec"]
                self.psm_dict["z"]=attributes["assumed_charge"]

                spectrum_native_id = "N/A"
                if "spectrumNativeID" in attributes:
                    self.psm_dict["spectrumNativeID"] = attributes["spectrumNativeID"]
                else:
                    self.psm_dict["spectrumNativeID"] = spectrum_native_id

                self.psm_dict["precursor_neutral_mass"] = attributes["precursor_neutral_mass"]
            if tag == "search_hit":
                self.psm_dict["num_proteins"] = attributes["num_tot_proteins"]
                self.psm_dict["protein"] = attributes["protein"]
                self.psm_dict["peptide"] = attributes["peptide"]
                self.psm_dict["matched_ions"] = attributes["num_matched_ions"]
                self.psm_dict["total_ions"] = attributes["tot_num_ions"]
                self.psm_dict["calc_neutral_mass"] = attributes["calc_neutral_pep_mass"]
                self.psm_dict["peptide_prev_aa"] = attributes["peptide_prev_aa"]
                self.psm_dict["peptide_next_aa"] = attributes["peptide_next_aa"]
            if tag == "modification_info":
                self.psm_dict["modification_info"] = attributes["modified_peptide"]
                #print("2",repr(self.psm_dict))
            if tag == "search_score":
                if attributes["name"] == "expect":
                    self.psm_dict["expect"] = attributes["value"]
            if tag == "alternative_protein":
                alternative_prots = []
                if "protein_alt" in self.psm_dict:
                    alternative_prots = self.psm_dict["protein_alt"]
                alternative_prots.append(attributes["protein"])
                self.psm_dict["protein_alt"] = alternative_prots
                #print("alt",attributes["protein"])
            if tag == "peptideprophet_result":
                self.psm_dict["pp_prob"] = attributes["probability"]
            if tag == "interprophet_result":
                self.psm_dict["ip_prob"] = attributes["probability"]
            if tag == "mod_aminoacid_mass" and "variable" in attributes:
                mod_aminoacid_masses = {}
                #print("1",repr(self.psm_dict))
                if "mod_aminoacid_mass" in self.psm_dict:
                    mod_aminoacid_masses = self.psm_dict["mod_aminoacid_mass"]

                mod_aminoacid_masses[attributes["position"]] = attributes["variable"]
                self.psm_dict["mod_aminoacid_mass"] = mod_aminoacid_masses




        def endElement(self, name):
            if name == "spectrum_query":
                if float(self.psm_dict["pp_prob"]) > PEPTIDE_PROPHET_THRESHOLD:

                    ppm_error = (float(self.psm_dict["precursor_neutral_mass"]) - float(self.psm_dict["calc_neutral_mass"]))*1000000/float(self.psm_dict["precursor_neutral_mass"])
                    da_error = float(self.psm_dict["precursor_neutral_mass"]) - float(self.psm_dict["calc_neutral_mass"])

                    output.write(self.psm_dict["spectrum"] + sep +
                                 self.psm_dict["spectrumNativeID"] + sep +
                                 self.psm_dict["precursor_neutral_mass"] + sep)
                    output.write(self.psm_dict["calc_neutral_mass"] + sep)
                    output.write(str(ppm_error) + sep)
                    output.write(str(da_error ) + sep)
                    output.write(self.psm_dict["z"] + sep)
                    output.write(self.psm_dict["rt"] + sep)
                    output.write(self.psm_dict["peptide"] + sep)

                    if "modification_info" in self.psm_dict:
                        output.write(self.psm_dict["modification_info"] + sep)
                    else:
                        output.write(self.psm_dict["peptide"] + sep)

                    if "mod_aminoacid_mass" in self.psm_dict:
                        mod_aminoacid_mass = self.psm_dict["mod_aminoacid_mass"]
                        mods = []
                        for pos in mod_aminoacid_mass:
                            mass = mod_aminoacid_mass[pos]
                            mods.append(mass + "@" + pos)
                        output.write(";".join(mods) + sep)
                    else:
                        output.write("" + sep)

                    output.write(self.psm_dict["peptide_prev_aa"] + sep +
                                 self.psm_dict["peptide_next_aa"] + sep)

                    output.write(self.psm_dict["pp_prob"] + sep)
                    ip_prob = ""
                    if "ip_prob" in self.psm_dict:
                        ip_prob = self.psm_dict["ip_prob"]
                    output.write(ip_prob  + sep)

                    output.write(self.psm_dict["matched_ions"] + sep +
                                 self.psm_dict["total_ions"] + sep)
                    output.write(self.psm_dict["protein"] + sep)
                    if "protein_alt" in self.psm_dict:
                        output.write(";".join(self.psm_dict["protein_alt"]))
                    #else:
                    #    output.write(sep)
                    output.write("\n")

    filename = os.path.basename(input)
    outfolder = os.path.dirname(input)
    #print(basename)

    if outfolder:
        outfolder = outfolder + "/"

    output_file = outfolder + file_prefix + "_" + filename[:-4]+".tsv"    #get rid of .pep.xml prefix, could mangle other things
    output = open(output_file, "w")
    output.write("spectrum" + sep + "native_id" + sep + "pre_neutral_mass" + sep + "calc_neutral_mass" + sep)
    output.write("ppm_error" + sep + "da_error" + sep)
    output.write("z" + sep + "rt" + sep + "peptide" + sep + "mod_peptide" + sep + "mod_masses" + sep + "pre" + sep + "post" + sep + "pp_prob" + sep + "ip_prob" + sep)
    output.write("matched_ions" + sep + "total_ions" + sep + "protein" + sep + "alt_protein" + "\n")


    # create an XMLReader
    parser = xml.sax.make_parser()
    Handler = TPPHandler()
    parser.setContentHandler(Handler)
    parser.parse(input)
    output.close()
    print("Finished. Output written to",output_file)



def testing():
    #input_XML = "D:/temp/PXD028712/interact-ipro.pep.xml"
    input_XML = "C:/Users/jonesar/Dropbox/DocStore/GrantApps/2025_BBSRC_PTMs_Var/UniProt_download/PXD010154/comet_materials/interact-ipro.pep.xml"
    #outfolder = "D:/temp/PXD028712/"
    OUT_PREFIX = "PXD010154"
    convert(input_XML,OUT_PREFIX)


if len(sys.argv)!= 3:
    print("Exit - expected usage\npython ",  sys.argv[0]," interact-ipro.pep.xml filenameprefix\n"
          "filenameprefix would usually be the PXD code or similar ")
else:
    convert(sys.argv[1],  sys.argv[2])

