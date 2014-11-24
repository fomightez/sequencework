#! /usr/bin/env python
 
# LookUpTaxonFA.py   *******by Wayne Decatur*********
# USES Python 2.7 because latest version on Python Anywhere with Biopython installed
# version 1.1
 
# FA version 1.1---> uses a list of FASTA entries from GenBank or the PDB to look up
# and assign Taxon (the Class specifically). Adds Taxon to end of description line
# for each FASTA entries. It also outputs to the user a report on the Taxons assigned.
# The idea is the FASTA entries have been already selected as to represent only
# a set of orthologs.
# Version 1.1 adds handling duplicates by not saving fasta entries if already
# have that organism's scientific name in list, as well as removing entries
# with duplicate accessions from output.  Also adds better tracking of
# taxon lookups to zero in on problem entries, like 'Chilabothrus striatus' faster.
#
# example of input and output:
# before:
# >gi|45383616|ref|NP_989588.1| follitropin subunit beta [Gallus gallus]
# MKTLNCYVLLFCWKAICCYSCELTNITIAVEREECELCITVNATWCSGYCFTRDPVYKYPPVSSVQQICT
# FKEVVYETVKIPGCGDHPESFYSYPVATECHCETCDTDSTDCTVRGLGPSYCSFSHNGSNQ
#
# after:
# >gi|45383616|ref|NP_989588.1| follitropin subunit beta [Gallus gallus]|Aves
# MKTLNCYVLLFCWKAICCYSCELTNITIAVEREECELCITVNATWCSGYCFTRDPVYKYPPVSSVQQICT
# FKEVVYETVKIPGCGDHPESFYSYPVATECHCETCDTDSTDCTVRGLGPSYCSFSHNGSNQ
#
# Might think could make first part faster and more direct by using the scientific name
# that is generally part of fasta entries instead of looking that up again using
# accession but 1) You'd be wrong because a lot of fasta description lines don't have
# that even though the majority do, and 2) foundation program (BIOALIGNAMER, aka namer)
# I had already went through looking
# this up as start step to getting records and luckily I left it not realizing I needed it.
# That program had been done by me because while I might remove a scientific name
# in interest of shorter description, I'd keep the accession info or would be able
# to get. Plus often fasta entries from other places don't have mcuh info and would
# be easier to build in simple accession number.
# I created it this way because I already had a program (BIOALIGNAMER, aka namer) that
# mined out some different taxonomy data using standard FASTA entries. Plus when I set
# to write this I wanted to BIOPYTHON to access NCBI taxonomy records as the method
# is optimized and because at the time on PythonAnywhere I couldn't access the GBIF
# http://data.gbif.org/ws/rest/taxon the Global
# BIODIVERSITY INFORMATION FACILITY's Taxon Web Service because it was not on the whitelist.
# *** UPDATE: AFTER I FOUND THE API for Global BIODIVERSITY INFORMATION
# FACILITY's Taxon Web Service I ASKED PYTHONANYWHERE TO ADD THIS SERVICE
# TO THE WHITELIST OF SITES FREE ACCOUNTS CAN ACCESS AND THEY DID!!! SO
# NOW I WILL INTEGRATE THIS CODE BACK INTO LookUpTaxonFA.py AND THE RELATED
# PROGRAMS AND CAN RUN THE SEARCH ALL AS ONE PROGRAM. AN IMPORTANT THING TO
# CONSIDER IS THAT ALTHOUGH THIS SYSTEM DEVELOPED HERE TO ACCESS THE TAXON DATA
# FROM the GBIF IS ALMOST FOOL-PROOF NOW, IT IS STILL BEST TO USE IT AS A
# REDUNDANT FALL BACK TO THE PRIMARY METHODS VIA BIOPYTHON SYSTEM ACCESSING NCBI.
# BIOPYTHON HAS BUILT IN ABILITY TO BE OPTIMIZED FOR HANDLING LARGE GROUPS
# AND I BUILD MY CODE TAKING ADVANTAGE OF THOSE ABILITIES. HOWEVER, IN THE CASE
# OF THE GBIF SYSTEM AND CODE I DEVELOPED FOR ACCESSING GBIF, THESE ARE NOT OPTIMIZED
# FOR LARGE NUMBERS OF LOOKUPS AND IN FACT MY METHOD IS BUILT TO BE THOROUGH
# AT THE COST OF BEING INEFFICIENT. SO USE BIOPYTHON-RELATED CODE TO GET MOST
# OF THE INFORMATION AND CLEAN UP THE TROUBLESOME ONES WITH THE CODE HERE TO
# ACCESS DATA AT THE GBIF BECAUSE THERE WILL BE ISSUES. SEVERAL ISSUES I HAD ALREADY
# CHARACTERIZED AND BECAUSE I DIDN"T EVEN TRY THAT MANY IN THE TESTING I DID TO
# DEVELOP MY TAXON FINDING CODE, I ASSUME THERE WILL BE MORE. THOSE I CHARACTERIZED
# AT THIS TIME ARE 1) THE NCBI DOESN'T LIST REPTILIA CLASS AT ALL FOR ANY REPTILES
# AND 2) DOES NOT LIST A CLASS FOR LAMPREY. (NCBI has Petromyzontiformes as ORDER
# WHICH I HAVE SEEN OTHERS USE AS THE TAXON FOR THIS ORGANIMS WHEN FOR OTHER
# ORGANIMS THE CLASS IS USED).
# ***END OF UPDATE***
#
#
#
#
#
# Related programs (besides namer v1.0.py, aka BIOALIGNAMER):
# FAComD version 1.0. (with FAComD in name) will --> start with list of MAFFT aligned modified
# fasta entries (modified to be aligned and have Common names via BIOALIGNAMER EARLIER). This will extract
# the numbers of the class Stacia needs for her ICCE slides. If I submit each group
# of sets of protein as a hand-copied block using the orresponding NEXUSbyPY file to
# tell me the boundaries of the six sets.
# While it would be nice to eventually build it so this program version
# accesses both files itself given the two names. That way it also opens the corresponding NEXUSbyPY
# file and mines the condensed, uniform protein designations which are appended at the end of the description
# line and that way it could compile a list automatically and the order is the same at
# that point so as long as I keep them in order they would match up (NOTE GIVEN I DIVIDE UP
# ENTRIES INTO GENBANK OR PDB SOURCE AS READ IN RIGHT NOW, THIS ADDITION OF ANOTHER FILE
# WOULD TAKE SOME CAREFUL TIMING AS DEALING WITH BOTH AND SO PROBABLY BEST TO MAKE A
# A STRING OF ALL LINES FOR EACH FILE AND THEN WORK THROUGH EACH IN PARaLLLEL.)
# Best to worry about PDB sourced files later because right now I don't have them
# and need the list for Stacia's slide. It would have been easier
# to already have that if I had added the condensed protein name to the description line
# at an earlier point in the workfow but I can work around not having the
# uniform protein designation for now without needing to do everything automatically a program by just breaking up
# the list by the corresponding sets of the glycoprotein.
#
#
# LookUpTaxonFORTABLE will --> start with edited list of description lines from fasta entries
# since I have this for forming a table with the necessary information
# to make a table like Uchida et al 2010 supplemental figure S1 (PMID: 20733079
# ) for my phylogenetic study of the glycoprotein hormones to go in Taka's
# manuscript. This program will add look up the class from NCBI's site for
# organisms extracted from what were previously the description lines for the
#FASTA entries before I altered them for the phylogenetic analysis and then
# altered them farther using REGEX as preliminary work towards this table.
# Because I need to tally the various vertebrate taxon for a slide in Stacia's
# ICCE 2013 presentation, I am going to add that as I go along.
 
 
##################################
#  USER SET VALUES               #
#REQUIRED ARGUMENTS BELOW
PathtoFolderWithFile = ""                                  #<---is for when running on PythonAnywhere.com
#PathtoFolderWithFile = "/Users/sowerlabwayne/Downloads/"  #<---is for when running local on my computer
 
#FileName= "tshbs.faa"
FileName = "testTAXON.faa"
#FileName = "testMMtable.txt"  #<---will be for version 1.2 without FA in name
 
UserEmail = "YOUR_EMAIL_HERE@example.com"
 
#OPTIONAL BELOW - STRUCTURED FOR ULTIMATELY USER ENTERED DATA SO ONLY SOME WILL GET ELEVATED TO FLAGS
CommonNamePreferred = False  # MOOT HERE
 
ProteinEntries = True  #set default TRUE
 
BothNamesAddedALL = True #OnlyComesUp AS USER-ENTERED DATA If CommonNamePreferred = True#Does User want Scientific name added in addition to common?
#SET DEFAULT FALSE SINCE MOST HAVE ALREADY IN FASTA #This value will enable adding both common and scientific names
#Tell user in most cases this is redundant as FASTA entry usually has Scientific name in header
 
BothNamesAddedJustPDB = True #ONLY TO COME UP AS USER-ENTERED DATA IF CommonNamePreferred = True, ProteinEntries = True and BothNamesAddedALL = False #SET DEFAULT TRUE
#This value will allow adding
 
TypicalOrder = True #ONLY TO COME UP AS USER-ENTERED DATA If CommonNamePreferred = True AND IF EITHER BothNamesAddedALL or BothNamesAddedJustPDB = TRUE
#Does user want typical order where common name listed first?  #Remind user this which is usually what gets used as designation in alignment and downstream work flows
#SET DEFAULT = TRUE
 
ModifyCaseOfCommonName = False #ONLY TO COME UP AS USER-ENTERED DATA IF CommonNamePreferred = True, #SET DEFAULT as false
#only use above to ask about additonal options - thus it should not come up in tests later
 
CapitalizeAllLettersOfCommonName = False #ONLY TO COME UP AS USER-ENTERED DATA IF CommonNamePreferred= True, #SET DEFAULT as false
 
CapitalizeFirstLetterOfCommonName = False #ONLY TO COME UP AS USER-ENTERED DATA IF CommonNamePreferred = True and CapitalizeAllLettersOfCommonName = False,
#SET DEFAULT as false
 
CapitalizeSubsequentLettersOfCommonName = False #ONLY TO COME UP AS USER-ENTERED DATA IF CommonNamePreferred = True
# and CapitalizeAllLettersOfCommonName = False, and CapitalizeFirstLetterOfCommonName = True    #SET DEFAULT as false
 
RemoveSpacesCommonName = False #ONLY TO COME UP AS USER-ENTERED DATA IF CommonNamePreferred#Does user want all whitespace removed from Common name so added
#as a single text value? #Tell users this sometimes helps with dealing with programs in your workflow that drop all but first word from consideration in designation
#SET DEFAULT as false BUT I WONDER IF SHOULD BE TRUE??? Maybe look into behavior of some of them to decide better.
 
RemoveSpacesScientifName = False #ONLY TO COME UP AS USER-ENTERED DATA IF CommonNamePreferred = false OR if BOTH (CommonNamePreferred AND BothNamesAddedALL= TRUE )
# OR if BOTH (CommonNamePreferred AND BothNamesAddedJustPDB= TRUE ) #Does user want all whitespace removed from Common name so added
#as a single text value? #Tell users this sometimes helps with dealing with programs in your workflow that drop all but first word from consideration in designation
#SET DEFAULT as false BUT I WONDER IF SHOULD BE TRUE??? Maybe look into behavior of some of them to decide better.
 
DashSpacesCommonName = False #ONLY TO COME UP AS USER-ENTERED DATA IF CommonNamePreferred= TRUE and RemoveSpacesCommonName = FALSE#Does user want all whitespace
# replaced by dashes in Common name so added as a single text value? #Tell users this sometimes helps with dealing with programs in your workflow that drop all
# but first word from consideration in designation #SET DEFAULT as false
 
DashSpacesScientifName = False #ONLY TO COME UP AS USER-ENTERED DATA IF CommonNamePreferred = false OR if BOTH (CommonNamePreferred AND BothNamesAddedALL= TRUE )
# OR if BOTH (CommonNamePreferred AND BothNamesAddedJustPDB= TRUE )#Does user want all whitespace replaced by dashes in Common name so added as a single text value?
#Tell users this sometimes helps with dealing with programs in your workflow that drop all but first word from consideration in designation #SET DEFAULT as false
##################################
 
 
 
 
#MAKE ALL ABOVE DOABLE WITH ARGUMENT FLAGS FIRST AND THEN ONLY ASK ABOUT REQUIRED IF NOT GIVEN, AND LATER ADD USER-ENTERED DATA MODE IF NO ARGUMENTS AND EVEN
#THEN ALLOW OPTION TO RUN WITH DEFAULTS SUBSEQUENTLY TO LET USER TEST
 
#BE SURE TO TELL USER IT WILL TAKE AT LEAST 5 SECONDS PER ENTRY TO WORK THROUGH
 
 
 
 
 
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###
#PART 1 - Getting user input on file and parameters and Reading in file designated by user
#
# for now this will be defined file
#TODO later ask user for file name to read in
#TODO To be added later is to also ask if want scientific name or common name added after greath-than symbol at start of fasta entry format for use in aligning DEFAULT IS COMMON NAME
 
import os
import sys
from Bio import Entrez
import urllib
import re
from collections import Counter
import time   #used for sleep
#import gzip
 
 
rootpath = PathtoFolderWithFile
fastalistfile = open(rootpath + FileName, "r")
NumberofFastaEntries = 0 #initiate with zero as value of Number of FASTA entries indentified
modifiedfastalistfilestring = "" #initiate a string that will become the modified file in the end
ListofAccessions =[] #initiate a list for the accessions
ListofPDBIds = [] #initiate a list for the PDBids,if any
DictionaryofIdsAndOldAndNewData = {} #initiate a dictionary to store Accession Id (or PDBid) as keys with value as a list containing placeholder and the DATA found to add
 
def SuggestReportBug():
    sys.stderr.write("\n         Program may work for you; however, you should be on the lookout for aberrations as the conditions triggering this WARNING should not occur. ")
    sys.stderr.write("\n         Please file a bug report concerning this warning with Wayne (send to wdecatur [AT] 'ye olde' yahoo.com ) and")
    sys.stderr.write("\n         please paste in or attach your input file unless it is mammoth or super secret. Any output to console you can include would help too.")
    sys.stderr.write("\n         \n")
 
def fetch_pdbheader(id):    #based on http://boscoh.com/protein/fetching-pdb-files-remotely-in-pure-python-code and http://www.pdb.org/pdb/static.do?p=download/http/index.html
    url = 'http://www.rcsb.org/pdb/files/%s.pdb?headerOnly=YES' % id
    return urllib.urlopen(url).read()
 
 
 
def fetch_keypage(akey):
    #first pause so don't cause site to be slammed when I have a lot to get
    time.sleep(4) #so don't hit the site too fast with a lot of calls
    #get record for gibifkey from GBIF taxon web service
    taxonurl = 'http://data.gbif.org/ws/rest/taxon/get?key=%s&stylesheet=' % akey
    #print taxonurl
    return urllib.urlopen(taxonurl).readlines()
 
 
def Fetch_gbifTaxon(taxonPageHandlerLinesForFunction,scinameinfo,keyinfo):
    extractedtaxon = "NotFound"  #intialize like this so will trigger alert below if never obtain
    #Now read through lines to find signal '<tn:rankString>class</tn:rankString>' and then want previous line
    #Then extract that Taxon from the previous line
    # TURNS OUT THE INSTANCES OF WEB PAGE SCRAPING NEEDED HERE WERE PRETTY BASIC BUT LOOKS LIKE IN FUTURE MAYBE WORTH LOOKING INTO BEAUTIFULSOUP http://stackoverflow.com/questions/7866128/python-split-without-removing-the-delimiter?lq=1 mentions it --> (v3 info because that is what Enthought allows ----> http://www.crummy.com/software/BeautifulSoup/bs3/documentation.html Pythonanywhere has 3 and 4)
    PreviousLine = "placeholder" #initialize with a placeholder in case something weird happens so it won't throw an error if needs it
    LineWithTaxa = "NotFound" #initialize this way so can trigger an alert if never found
    MatchToSearchFor ='<tn:rankString>class</tn:rankString>'
    for line in taxonPageHandlerLinesForFunction:
        #print line
        if MatchToSearchFor in line:
            #get last line and extract GBIF species number
            LineWithTaxa = PreviousLine
            break #because no reason to keep cycling through lines of HTML if have it , see http://www.tutorialspoint.com/python/python_loop_control.htm ('continue' statement there looks interesting for situations and I could use 'pass' a lot in writing when need to leave a stub)
        PreviousLine = line
    if LineWithTaxa == "NotFound":  #giving feedback of what steps seemed to fail, if any
        #report warning alert because should be able to determine taxa of scientific name correct
        sys.stderr.write("\n\n*******************************************************************************")
        sys.stderr.write("\n*******************************ALERT***ALERT***********************************")
        sys.stderr.write("\n******      Unable to find the taxon in the GBIF record returned    *********** ")
        sys.stderr.write("\n*******     for "+ scinameinfo +" using gbifkey of "+keyinfo+" .    *******")
        sys.stderr.write("\n*******************************************************************************")
        sys.stderr.write("\n*******************************************************************************\n\n")
    #Extract the taxon from the line that has it = LineWithTaxa
    #taxaline should look something like '<tn:nameComplete>Reptilia</tn:nameComplete>'
    #Now need to get the actual class out of between <tn:nameComplete>...</tn:nameComplete>' tags by deleting the tags
    extractedtaxon = LineWithTaxa.strip().replace('<tn:nameComplete>','').replace('</tn:nameComplete>','')
    if extractedtaxon == "NotFound":
        #report warning alert because should be able to determine taxa of scientific name correct
        sys.stderr.write("\n\n*******************************************************************************")
        sys.stderr.write("\n*******************************ALERT***ALERT***********************************")
        sys.stderr.write("\n*********        Unable to extract taxon from the GBIF record of      *********")
        sys.stderr.write("\n*****         key "+keyinfo+" for "+scinameinfo+".        *****")
        sys.stderr.write("\n*******************************************************************************")
        sys.stderr.write("\n*******************************************************************************\n\n")
    return extractedtaxon
 
 
 
def fetch_taxonfromgbif(sciname):
    gbifkey = "NotFound"  #intialize like this so will trigger alert if never obtain
    ListofDeadEndgbifKeys = []
    gbifTaxon = "NotFound"  #intialize like this so will trigger alert if never works
    #because it will be part of a web address, need to %xx escape the scientific name for URL encoding. According to http://stackoverflow.com/questions/8905864/url-encoding-in-python it is easy.
    escapedSciName = urllib.quote(sciname)
    #print escapedSciName
    #WRONGurl = 'http://data.gbif.org/search/%s' % escapedSciName
    #DOESN"T WORK BECUASE JAVASCRIPT RENDERS NORMAL PAGES LIKE AT --->http://data.gbif.org/search/Petromyzon%20marinus
 
    # I has searched the GBIF site myself at first trying to find this but finally after realizing data was hidden I searched 'python gbif' which lead to
    # https://github.com/matagus/python-gbif which lead to the API at  ----> http://data.gbif.org/ws/rest/taxon
    '''
    parameter_list is a set of URL-encoded key-value pairs including any of
            the following keys:
 
          scientificname - return only records where the scientific name matches
              that supplied.
    '''
 
    '''
    Because about stylesheet "Supplying an empty string as the value for this parameter
              disables the setting of a stylesheet by the data portal." HOORAY!!!
    NOTE THOUGH THE PAGE IN A BROWSER WILL LOOK BLANK!!! VIEW SOURCE TO SEE TEXT!!!
    EXAMPLE
    http://data.gbif.org/ws/rest/taxon/list?scientificname=Crocodylus%20acutus&stylesheet=
 
    THEN
    Get species number from finding text like this, where 114987548 is the key:
        <tc:TaxonConcept gbifKey="114987548" status="unconfirmed"
 
    Alright once I have gbif key I should grab page with
    EXAMPLE
    http://data.gbif.org/ws/rest/taxon/get?key=114987548&stylesheet=
 
    And then locate on that page this information:
 
    <tn:nameComplete>Reptilia</tn:nameComplete>
         <tn:rankString>class</tn:rankString>
    '''
    url = 'http://data.gbif.org/ws/rest/taxon/list?scientificname=%s&stylesheet=' % escapedSciName
    #print url
    gbifListPageHandler = urllib.urlopen(url).read()
    #print gbifListPageHandler
    IndentifyingCodeToSearchFor = '<tc:TaxonConcept gbifKey="'
    #so can split the page using that match xml code
    ChoppedPage = gbifListPageHandler.split(IndentifyingCodeToSearchFor)
    #Now to get just the gbifKey, split the second entry in the ChoppedPage list
    # where quotes is and that should leave the key as first entry
    # TURNS OUT THE INSTANCES OF WEB PAGE SCRAPING NEEDED HERE WERE PRETTY BASIC BUT LOOKS LIKE IN FUTURE MAYBE WORTH LOOKING INTO BEAUTIFULSOUP http://stackoverflow.com/questions/7866128/python-split-without-removing-the-delimiter?lq=1 mentions it --> (v3 info because that is what Enthought allows ----> http://www.crummy.com/software/BeautifulSoup/bs3/documentation.html Pythonanywhere has 3 and 4)
    indexthatshouldbegbif = 1 #index to check first should be 1 but I am going to make a loop below for when doesn't work first time so going to set as a variable
    gbifkey=(ChoppedPage[indexthatshouldbegbif].split("\""))[0]  #backslash is to escape quote
    #print gbifkey
 
    #give feedback of what steps seemed to fail, if any
    if gbifkey  == "NotFound":
        #report warning alert because should be able to determine taxa of scientific name correct
        sys.stderr.write("\n\n*******************************************************************************")
        sys.stderr.write("\n*******************************ALERT***ALERT***********************************\n")
        sys.stderr.write("\n Unable to extract gbif key number at the GBIF site ")
        sys.stderr.write("for "+ sciname +".")
        sys.stderr.write("\n*******************************************************************************")
        sys.stderr.write("\n*******************************************************************************\n\n")
        return gbifTaxon   #no sense going farther since failed to get key from what was returned. Since using same method below for getting a different key, it won't help trying that again.
 
    #Now use this key to obtain the record that should have the Taxon (specifically class) on it
    if gbifkey != "NotFound":   #because doesn't make sense to try for page from GBIF taxon web service if the key wasn't found
        taxonPageHandler = fetch_keypage(gbifkey)
        #print taxonPageHandler
        #Cycle through data returned and Extract the taxon from the line that has it
        gbifTaxon = Fetch_gbifTaxon(taxonPageHandler,sciname,gbifkey)
        #print gbifTaxon
        if gbifTaxon == "NotFound":
            ListofDeadEndgbifKeys.append(gbifkey)
        else:
            sys.stderr.write("....key worked....taxon obtained... ")
            return gbifTaxon
 
 
    #Try again with a different key from that page until get something
    ExtraTriesPossible= len(ChoppedPage)/2 - 2  #Can only try if there are more keys from the record (which should be half the number of elements in the list because of the way it was chopped); the subtraction is to deal with the first key already used
    #print ExtraTriesPossible
    for x in xrange(ExtraTriesPossible):
        gbifkey = "NotFound"  #reintialize
        sys.stderr.write("... key failed....trying with a different key.... ")
        indexthatshouldbegbif += 2   #the next key should be 2 indices from previous
        gbifkey=(ChoppedPage[indexthatshouldbegbif].split("\""))[0]  #backslash is to escape quote
        #if that key is in the deadkeys list, we need another before trying
        if gbifkey in ListofDeadEndgbifKeys:
            continue     # that continue command did come in handy, see where I mentioned seeing this in http://www.tutorialspoint.com/python/python_loop_control.htm ; I was going to use a while loop to get another key until I had one not in list but that would potentially become an endless loop if all keys already in there, this is much better for staying within possible numbers
        #print gbifkey
 
        #Now use this key to obtain the record that should have the Taxon (specifically class) on it
        if gbifkey != "NotFound":   #because doesn't make sense to try for page from GBIF taxon web service if the key wasn't found
            taxonPageHandler = fetch_keypage(gbifkey)
            #print taxonPageHandler
            #Cycle through data returned and Extract the taxon from the line that has it
            gbifTaxon = Fetch_gbifTaxon(taxonPageHandler,sciname,gbifkey)
            #print gbifTaxon
            if gbifTaxon == "NotFound":
                ListofDeadEndgbifKeys.append(gbifkey)
            else:
                sys.stderr.write("....key worked....taxon obtained... ")
                return gbifTaxon
 
    return gbifTaxon  #this should never get used unless just cannot find SO KEEP HERE FOR THOSE POSSIBLE CASES
 
 
 
 
def FixCase(Genusspecies):
    # to assure case is correct first lower all and then make first upper (could have uppered all first and the did inverse)
    Genusspecies = Genusspecies.lower()
    return Genusspecies [0].upper()+Genusspecies [1:]
 
def GetTaxon(thesciname):
    #First try NCBI via BioPython
    Taxon = "NotFound" #initialize with this so we know if found anything
    if thesciname == DuplicatesIndicatorTag:
        Taxon = DuplicatesIndicatorTag
        return Taxon
    sys.stderr.write(".")
    if ('(' in thesciname) or (')' in thesciname):  #found problem with 'Xenopus (Silurana) tropicalis' not returning any
        # taxonomy entries and so thought might work if removed parentheses and it does, actuallly works better with 'Silurana tropicalis'
        # that works for http://data.gbif.org/  which may be important. I AM SWITCHING TO MAKING IT LAST TWO WORDS WITHOUT PARENTHESES
        ScientificNameForSearch=thesciname.replace('(','').replace(')','')
        words = ScientificNameForSearch.split() # so I can count and deal with words
        if len(words) >2:
            ScientificNameForSearch = ' '.join(words[-2:])  #takes just the last two words of the list
            ScientificNameForSearch = FixCase(ScientificNameForSearch)  # to insure the case is correct now
        #print ScientificNameForSearch
    else:
        ScientificNameForSearch=thesciname  #will be used for search. Most times no need to be modified.
    handle = Entrez.esearch(db="taxonomy", term=ScientificNameForSearch, usehistory="y")
    esearchrecord = Entrez.read(handle)
    webenv2 = esearchrecord["WebEnv"]
    query_key2 = esearchrecord["QueryKey"]
    #print esearchrecord
    #print esearchrecord["Count"]
    if int(esearchrecord["Count"])>1:
         sys.stderr.write("\nWARNING: More than one taxonomy record came up for "+  thesciname +"; however, the first one was used.")
         sys.stderr.write("\n         You should check the selected name and hand edit output at end, if needed.")
    if len(esearchrecord["IdList"]) ==0:
        sys.stderr.write("\nWARNING: No taxonomy record came up for "+  thesciname +".")
        sys.stderr.write("\n         You should check the selected name and hand edit output at end, if needed.")
    else:
        sys.stderr.write("..communicating with NCBI...looking up taxon of "+  thesciname +".")
        fetch_handle = Entrez.efetch(db="taxonomy", id=esearchrecord["IdList"][0], webenv=webenv2, query_key=query_key2)
        sys.stderr.write(".")
        #taxonomydatarecords = Entrez.read(fetch_handle)
        taxonomydata = fetch_handle.read()
        #print taxonomydata
        #I was actually successful in getting the BIOPYTHON Entrez parser to parse this as above but I found the resulting hierarchy too nested, i.e.,
        # it was a list of dictionaries that was to make heads or tails of easily and
        #this approach below using minidom that had researched earlier seemed more direct
        # however this example is much messier than the case with finding commmon name in Bioalignamer
        #but it workS for now. The clunky part is that I go through the LineageEx lines to get the line
        # before the '<Rank>class</Rank>' line. There should be a simple way to get that.
        from xml.dom.minidom import parseString
        dom = parseString(taxonomydata)
        if dom.getElementsByTagName('LineageEx').length > 0:     #<---------from
        # http://stackoverflow.com/questions/9091549/xml-dom-minidom-count-how-many-tags-of-a-specific-type-exist
        # and see http://bytes.com/topic/python/answers/699070-checking-whether-xml-dom-node-already-exists
        #NEED THAT ABOVE LINE TO MAKE SURE NO ERROR THROWN IF THE TAG DOESN'T EXIST
            xmlTagLinesFIRST = dom.getElementsByTagName('LineageEx')[0].toxml()   #Gets first instance of LineageEx and all the lines between its closing tag plus closing
            xmlTagLines = str(xmlTagLinesFIRST) #<---found I had to convert this to string and split it to be able to check if "<Rank>class</Rank>" in it below
            #print "This is XMLTAGLINES"+xmlTagLines
            xmlTagLines= xmlTagLines.split()
            #Now that have XML lineage need to step through lines and grab the line before '<Rank>class</Rank>' line
            #because that will have the CLASS I want! If it is there
            lastline = "placeholder" #initialize so on off chance first line has something
            LineContainsClass = "NotFound" #initialize with this in case not there so when Taxon assigned below it gets the alert signal
            for theline in xmlTagLines:
                if "<Rank>class</Rank>" in theline:
                    LineContainsClass = lastline.strip ()
                    sys.stderr.write("..Taxon found..")
                    #Now need to get the actual class out of between <ScientificName>...</ScientificName> tags by deleting the tags
                    Taxon = LineContainsClass.replace('<ScientificName>','').replace('</ScientificName>','')
                    return Taxon    #because return breaks the loop (see http://www.youtube.com/watch?v=K96chCMCQ68), putting it here will step reading the rest of the lines once it has Taxon, optimizes time
                lastline = theline
 
        # At this point if Taxon still is "NotFound" then either 1) the NCBI entry seems odd
        # and doesn't have LineageEx XML or 2) I was not able to grab the class from
        # list as it most likely doesn't exist as is the case with Reptilia.
        ''' Additional notes on the Reptilia class problem:
        NCBI taxonomy server doesn't list that classs
        as it is divided among three suborders(?):
            Lepidosauria [in-part: Reptilia]
            Testudines [in-part: Reptilia]
            Crocodylidae [in-part: Reptilia]
 
        BUT THEY AREN'T EVEN DISPLAYED CONSISTENTLY IN THE NCBI RECORDS
        ,for example see Anolis carolinensis and Crocodylus acutus . The conventions
        aren't consistent for either relative what the Reptilia are divided into and
        so if that is going to be annoying I'd rather devise a more general fall
        back place to work with such problematic taxonomic groups because:
        1) if it problems came up with the NCBI taxonmy with one of the first 5 I tested it with then there is bound to be others;
        and 2) it seems there is alternative one at the http://data.gbif.org/ . Because http://data.gbif.org/search/Anolis%20carolinensis
        would work as does http://data.gbif.org/search/Crocodylus%20acutus and the clearly show the class taxon there!!
        However, I am going to keep the primary as NCBI because it should work for most and is the better to do this then potentially
        hitting the Global BIODIVERSITY INFORMATION FACILITY a lot. This will be similar
        to how I handled the fasta files from the PDB in the bioalignamer program I wrote.
        '''
        # Either way, next try the http://data.gbif.org/ws/rest/taxon the Global
        # BIODIVERSITY INFORMATION FACILITY's Taxon Web Service. alternative
        if Taxon == "NotFound":
            #LOOK UP AND GET THE CLASS FROM GBIF taxon web service API as described at http://data.gbif.org/ws/rest/taxon
            sys.stderr.write("\n...Trying gbif.org for taxon..")
            Taxon = fetch_taxonfromgbif(thesciname)
            return Taxon
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#*********************************************************************************
#*********************************************************************************
#*********************************************************************************
#********MAIN PART OF PROGRAM STARTS HERE*************** ABOVE ARE FUNCTIONS******
sys.stderr.write("Reading in your FASTA entries..")
for line in fastalistfile:
    line = line.strip ();#this way I have better control of ends ultimately becaue I know there isn't any unless I add
    if len(line) > 0:
        if line.startswith('>'): # only need those that start with the greater-than symbol as they have the accessions and need to be modified
            NumberofFastaEntries += 1
            sys.stderr.write(".")
            words = line.split()    # In case I need words as I go through
            InfoAndUIDs = line.split("|") #makes it easier to grab designation of type and accession number than words
            FormatOfFasta = (InfoAndUIDs[0])[1:]
            Accession = InfoAndUIDs[1] # Actually because it is easy to get and works I am using the GI number, seee http://www.ncbi.nlm.nih.gov/Sitemap/sequenceIDs.html
            FormatAndAccession = FormatOfFasta + '|' + Accession #removes the ">" from beginning of Format string and appends accession
            #Because Fasta files coming from PDB have a different arrangement, need to adjust and want to add to a different list
            if (Accession == 'PDBID'):  #Accession object not accession value because has format from RCSB output so need to fix
                FormatOfFasta = Accession
                Accession = (InfoAndUIDs[0])[1:5] #RCSB has PDB id code here in PDB FASTA
            #print InfoAndUIDs
            #print FormatAndAccession
            FormatAndAccession = FormatOfFasta + '|' + Accession  #Fix this now that set correct one # I DON'T THINK USED FURTHER
            Placeholder = 'PLACEHOLDER'+FormatAndAccession
            DuplicatesIndicatorTag = "DuplicateAlreadyInList"
            if (InfoAndUIDs[1] == 'PDBID'):
                if Accession not in ListofPDBIds:  #Do not bother adding if already there because it will save us time later by avoiding repeating lookups
                    ListofPDBIds = ListofPDBIds + [Accession]
                    line = '\n'+ line + Placeholder     #want placeholder at end so can tack on TAXON at end
                else: #if it is a duplicate ID
                    #if already seen tag so in addition to not adding to lookup list it will later be deleted when script handles removing the duplicates from the list of fasta entries
                    line = '\n'+ line + " "+ DuplicatesIndicatorTag+" " + Placeholder     #want placeholder at end so can tack on TAXON at end EVEN THOUGH THIS LINE SHOULD GET REMOVED AS A DUPLICATE ULTIMATELY
            else:
                if Accession not in ListofAccessions:  #Do not bother adding if already there because it will save us time later by avoiding repeating lookups
                    ListofAccessions = ListofAccessions + [Accession] #this will later be submitted Entrez as at http://www.bio-cloud.info/Biopython/en/ch8.html
                    line = '\n'+ line + Placeholder     #want placeholder at end so can tack on TAXON at end
                else: #if it is a duplicate ID
                    #if already seen tag so in addition to not adding to lookup list it will later be deleted when script handles removing the duplicates from the list of fasta entries
                    line = '\n'+ line + " "+ DuplicatesIndicatorTag+" " + Placeholder     #want placeholder at end so can tack on TAXON at end EVEN THOUGH THIS LINE SHOULD GET REMOVED AS A DUPLICATE ULTIMATELY
            DictionaryofIdsAndOldAndNewData [Accession] = [Placeholder,"NotLookedUp"]    #NotLookedUp as a marker of status
            #print Accession
        modifiedfastalistfilestring = modifiedfastalistfilestring + '\n' + line #sequence lines can just be passed to what will make of the new file
#print modifiedfastalistfilestring
#print ListofAccessions
#print ListofPDBIds
#print DictionaryofIdsAndOldAndNewData.values()[0][1]
 
#just a check that all is well when comparing Accession/PDB IDs parsed and the dictionary used for storing placeholder and name, OR notify user if kosher.
#SHOULD NOT TRIGGER.
if (len (ListofAccessions) + len(ListofPDBIds)) != len (DictionaryofIdsAndOldAndNewData):
    sys.stderr.write("\nWARNING: Conflict in total number of combined Accessions and PDBIds read in as compared to the total number stored in list keeping track of names.")
    SuggestReportBug()
 
#give user some stats
sys.stderr.write(" Concluded. \n"+ str(NumberofFastaEntries)+ " FASTA entries read in with " + str(len (ListofAccessions))+ " unique Accession numbers")
if len(ListofPDBIds)>0:
    sys.stderr.write(" and " + str(len (ListofPDBIds))+ " PDB IDs.")
else:
    sys.stderr.write(".")
if NumberofFastaEntries != ListofAccessions:
    sys.stderr.write("\nThe entries with duplicate accessions will be removed in the output.")
#done with user file read in
fastalistfile.close()
 
 
 
 
 
 
 
 
 
#get full entries to extract the Scientific names
#FIRST DOING ALL BUT PDBs
sys.stderr.write("\nObtaining standardized records related to entries:")
# Doing it this way because while many entries, such as GI: 18640742, have source organism scientific name (genus and species) at end of entry,
# many do not, such as GI:121313. Thus looking up will be more versatile.
#Also EVENTUALLY makeit so it tests if each entry is protein or nucleic acid THIS SHOULD BE ADDED AS ANOTHER VALUE TO THE DICTIONARY
#AND CAN BE DONE AS THE NEXT LINE after header  IS READ IN BECAUSE THE the accession will still be in memory (THIS WOULD HAVE A PROBLEM WITH FASTA FILES WITH DOUBLE HEADERS THOUGH BUT PROBABLY RARE!!)
 
#initiate lists and strings that get added to
ScientificNameList = []
 
#submit all as a block as that is the best way to let biopython and NCBI manage it according to http://www.bio-cloud.info/Biopython/en/ch8.html
if len(ListofAccessions)>0:
    Entrez.email = UserEmail
    if ProteinEntries:
        handle = Entrez.epost("protein", id=",".join(ListofAccessions))
    else:
        handle = Entrez.epost("nucleotide", id=",".join(ListofAccessions))
    ePostResult = Entrez.read(handle)
    webenv = ePostResult["WebEnv"]
    query_key = ePostResult["QueryKey"]
    #print ePostResult
    FullEntries = ""
    count = len(ListofAccessions)
    batch_size = 20
    for start in range(0,count,batch_size):
        end = min(count, start+batch_size)
        sys.stderr.write("...downloading records %i to %i..........................." % (start+1, end))
        if ProteinEntries:
            handle = Entrez.efetch(db="protein", rettype="gb",
                                     retstart=start, retmax=batch_size,
                                     webenv=webenv, query_key=query_key)
        else:
            handle = Entrez.efetch(db="nucleotide", rettype="gb",
                                     retstart=start, retmax=batch_size,
                                     webenv=webenv, query_key=query_key)
        FullEntries = FullEntries + handle.read()
 
    #SINCE DOWNLOADING AS A GROUP CANNOT JUST SEND EACH ONE TO A LIST SO FOUND A WAY TO BREAK BACK UP TO A LIST
    #by the way I probably could have just used '//' that is at the end of each GENBANK entry but didn't because was worried
    # since '/' used a lot in entries it might show up other than at end in entries and I didn't what to have to invoke regular expression to insure at start of a line AND I sould still have to deal with // at end.
    #In hindsight though the approach of using '//' might have been easier and probably my question if it can occur in middle might be answered if I looked into format of genpept records.
    FullEntriesList = FullEntries.split("LOCUS       ")
    del FullEntriesList[0] #removes first table element because 'LOCUS       ' separator used to indentify split occurs at very start of string so Python adds an extra empty element in front of it
    #Next line puts back "Locus" that got deleted by split, using re.split leaves in list but not in location so this approach ended up being easiest
    FullEntriesList = ["LOCUS       "+FullEntriesList[i] for i in range(len(FullEntriesList)) ]  #using list comprehension as described http://desk.stinkpot.org:8080/tricks/index.php/2007/10/extract-odd-or-even-elements-from-a-python-list/
 
 
 
    #NOW MINE FOR THE NAMES
    NumberofSciNameDuplicatesNoted = 0 #initialize Duplicate tracking valiue
    sys.stderr.write("\nMining scientific ")
    if CommonNamePreferred:
        sys.stderr.write("and common ")
    sys.stderr.write("names of source organisms..")
    for entry in FullEntriesList:
 
        #print entry[0:2000]
        ScientificName="" #Resetting so no bleed through during loops
        GenbankCommonNameData="" #Resetting so no bleed through during loops
        CommonNameData="" #Resetting so no bleed through during loops
 
        sys.stderr.write(".")
        ScientificNameCatcherList = re.split ('  (ORGANISM)',entry)
        ScientificName = ((ScientificNameCatcherList[ScientificNameCatcherList.index('ORGANISM')+1].splitlines())[0]).strip() #grabs
        # the next element in the OrganismCatcherList after 'ORGANISM' and splits it into lines and returns the first line,
        # plus strips off leading and trailing spaces I found showing up
        #AS AN AUTOMATED MEASURE TO DEAL WITH DUPLICATES, NEXT ONLY ACTUALLY ADD
        # NAME IF DON'T ALREADY HAVE THAT ONE REPRESENTED THERE ALREADY,
        # otherwise add a recognizable tag signifying that and I can deal with it
        # later
        if ScientificName not in ScientificNameList:
            ScientificNameList.append(ScientificName)
        else:
            ScientificNameList.append(DuplicatesIndicatorTag)
            #report warning alert because should be able to determine taxa of scientific name correctly
            sys.stderr.write("\n*******************************NOTE***NOTE*************************************\n")
            sys.stderr.write("\n ***  DUPLICATE.  "+ ScientificName +" already in list.   ***")
            sys.stderr.write("\n*******************************************************************************\n\n")
            NumberofSciNameDuplicatesNoted += 1
        #print ("the organism name is "+ ScientificName +".")
    #print ScientificNameList
 
 
 
 
#Now to deal with PDB FASTAs to get Scientific Names - can take shortcut of just getting headers of PDB entries since include
sys.stderr.write("\n")
PDBScientificNameList = [] #by putting here, later lines won't return error using it to make total list
if len(ListofPDBIds)>0:
    sys.stderr.write("...Obtaining record related to a PDB entry.....")
    #initiate lists and strings that get added to
    PDBPageList = []
    for PDBid in ListofPDBIds:
        if ListofPDBIds.index(PDBid) != 0:  #for every one after the first, want pause
            time.sleep(4) #so don't hit the site too fast with a lot of calls
        #sys.stderr.write(".")
        PDBScientificName="" #Resetting so no bleed through during loops
        PDBCommonName="" #Resetting so no bleed through during loops
        PDBhandler = fetch_pdbheader(PDBid)
        #print PDBhandler[0:2800]
        #PDBPageList.append(PDBhandler) #DECIDED I DIDN'T NEED #couldn't pass each to a list entry with Full entries from NCBI because read in in batches so would add more than one anyway
        if 'ORGANISM_SCIENTIFIC:' in PDBhandler:
            sys.stderr.write(".")
            PDBScientificNameCatcherList = re.split ('(ORGANISM_SCIENTIFIC:)',PDBhandler)
            #print PDBScientificNameCatcherList
            PDBScientificName = ((PDBScientificNameCatcherList[PDBScientificNameCatcherList.index('ORGANISM_SCIENTIFIC:')+1].splitlines())[0]).strip() #grabs
            # the next element in the PDBScientificNameCatcherList after 'ORGANISM_SCIENTIFIC:' and splits it into lines and returns the first line,
            # plus strips off leading and trailing spaces I found showing up
 
            PDBScientificName = PDBScientificName[0]+PDBScientificName[1:-1].lower()#since seems all caps from PDB, made all but first letter lowercase to be more like normal
            # and I added command to take only up to last character since was coming backwith semi-colon at end of name
            #print "\nPDBScientific NAme = "+ PDBScientificName +"\n"
        else:
            #report warning because should be there
            sys.stderr.write("\n\n*******************************************************************************")
            sys.stderr.write("\n*******************************ERROR***ERROR***********************************\n")
            sys.stderr.write("\nSorry for any problem this causes but....\n Unable to find scientific name for source \
            organism for")
            sys.stderr.write(" PDB entry "+ PDBid +".")
            sys.stderr.write("\n*******************************************************************************")
            sys.stderr.write("\n*******************************************************************************")
            PDBScientificName = 'NoneObtained'
        PDBScientificNameList.append(PDBScientificName)
 
 
 
 
 
#Now to go through the two lists and get the Taxon, specifically THE CLASS information
#For most, NCBI via BIOPYTHON will suffice but I saw problem with reptiles so I'll use
# http://data.gbif.org/ws/rest/taxon the Global BIODIVERSITY INFORMATION FACILITY to get ones
# that don't list a class entry in XML Taxonmoy records
# Do ScientificNameList first
# initialize list for taxon (class)
taxonlist = []
NotFoundTaxaList = []
TheCurrentTaxon = "NotFound" #initialize and insure it triggers and alert if it doesn't work
sys.stderr.write("\n")
sys.stderr.write("Looking up taxons for the entries..")
for Orgname in ScientificNameList:
    if (ScientificNameList.index(Orgname) != 0) and (Orgname != DuplicatesIndicatorTag):  #for every one after the first, want pause, but don't need if no scientific name being looked up because it is a duplicate
        time.sleep(4) #so don't hit the site too fast with a lot of calls
    #get Taxon for Orgname
    sys.stderr.write(".")
    TheCurrentTaxon = GetTaxon(Orgname)
    if TheCurrentTaxon == "NotFound":
        #report warning alert because should be able to determine taxa of scientific name correctly
        sys.stderr.write("\n\n*******************************************************************************")
        sys.stderr.write("\n*******************************ALERT***ALERT***********************************\n")
        sys.stderr.write("\n ***  Unable to find the taxon ")
        sys.stderr.write("for "+ Orgname +".   ***")
        sys.stderr.write("\n*******************************************************************************")
        sys.stderr.write("\n*******************************************************************************\n\n")
        NotFoundTaxaList.append(Orgname)
    taxonlist.append(TheCurrentTaxon)
 
# Next do PDBScientificNameList
PDBtaxonlist = []
TheCurrentTaxon = "NotFound" #reinitialize and insure it triggers and alert if it doesn't work
sys.stderr.write("Looking up taxons for your entries that came from the PDB..")
for PDBOrgname in PDBScientificNameList:
    if PDBScientificNameList.index(PDBOrgname) != 0:  #for every one after the first, want pause
        time.sleep(4) #so don't hit the site too fast with a lot of calls
    #get Taxon for PDBOrgname
    if PDBOrgname != 'NoneObtained':
        TheCurrentTaxon = GetTaxon(PDBOrgname)
    if TheCurrentTaxon == "NotFound":
        #report warning alert because should be able to determine taxa of scientific name correctly
        sys.stderr.write("\n\n*******************************************************************************")
        sys.stderr.write("\n*******************************ALERT***ALERT***********************************\n")
        sys.stderr.write("\n ***  Unable to find the taxon ")
        sys.stderr.write("for "+ PDBOrgname +".   ***")
        sys.stderr.write("\n*******************************************************************************")
        sys.stderr.write("\n*******************************************************************************\n\n")
        NotFoundTaxaList.append(PDBOrgname )
    PDBtaxonlist.append(TheCurrentTaxon)
 
 
 
 
 
 
 
#Tallies for Taxa
sys.stderr.write(" \n***********THE ANALYSIS IS COMPLETE************\n")
sys.stderr.write(" NotFoundTaxaList is \n"+ str(NotFoundTaxaList) +"\n")
CompleteTaxaList = taxonlist + PDBtaxonlist
OccurencesCompleteTaxaList = dict(Counter(CompleteTaxaList))
if DuplicatesIndicatorTag in OccurencesCompleteTaxaList:  #otherwise it throws an error if try to remove key on the next line and it doesn't exist
    del OccurencesCompleteTaxaList[DuplicatesIndicatorTag] #REMOVES DICTIONARY ARISING FROM THE USE OF DUPLICATES TAG , see http://stackoverflow.com/questions/5447494/best-way-to-remove-an-item-from-a-python-dictionary
sys.stderr.write(" Number of occurences for each Taxon is:\n"+ str(OccurencesCompleteTaxaList) +"\n")
 
 
 
 
 
 
 
#NOW TO MERGE RESULTS back into string of modified FASTA file
 
#first prepare element in position 2 of dictionary value so can know what got changed and what happened when
for TheKey in DictionaryofIdsAndOldAndNewData.keys():
    (DictionaryofIdsAndOldAndNewData[TheKey])[1] = "NotObtained"   #changing status after other step and before next makes easier to track progress later
 
#fix values in ScientificNameList and PDBScientificNameList to be put in according to optional parameters
#since want lists for final stats, make a new list with processed ones
ProcessedNamesList = ScientificNameList    #start with current list
PDBProcessedNamesList = PDBScientificNameList #start with current list
#process for ProcessedNamesList
if CapitalizeAllLettersOfCommonName:
    ProcessedNamesList = [Val.upper() for Val in ProcessedNamesList]
if CapitalizeFirstLetterOfCommonName and CapitalizeSubsequentLettersOfCommonName:
    ProcessedNamesList = [Val.title() for Val in ProcessedNamesList]
#process for PDBProcessedNamesList
if CapitalizeAllLettersOfCommonName and (len(ListofPDBIds)> 1):
    PDBProcessedNamesList = [Val.upper() for Val in PDBProcessedNamesList]
if CapitalizeFirstLetterOfCommonName and CapitalizeSubsequentLettersOfCommonName and (len(ListofPDBIds)> 1):
    PDBProcessedNamesList = [Val.title() for Val in PDBProcessedNamesList]
#print ProcessedNamesList
#print PDBProcessedNamesList
 
 
 
#NOW go through the accession list and use accession id as a key to select the dictionary entry to process and then change the text at list position 1 for so it is the TAXON
sys.stderr.write("....updating FASTA entries in output.....")
#next for statement will add to the dictionary entry what I'll want to swap into the current placeholder linked to each key id
for eachaccession in ListofAccessions:
    (DictionaryofIdsAndOldAndNewData[eachaccession])[1] = "|"+ taxonlist[ListofAccessions.index(eachaccession)]  # I added "|" in front of taxon for a separator
#do similar for PDBIDS
if len(ListofPDBIds)>0:
    sys.stderr.write("....updating PDB-related FASTA entries in output.....")
    for eachpdbid in ListofPDBIds:
        (DictionaryofIdsAndOldAndNewData[eachpdbid])[1] = "|"+ PDBtaxonlist[ListofPDBIds.index(eachpdbid)]  # I added "|" in front of taxon for a separator
#NEXT THE BIG ACTUAL REPLACEMENTS THAT MERGE TAXONS TO END OF EACH DESCRIPTION LINE
#now that each entry in dictionary has name value in its list set, cycle through all the keys and replace placeholder text (element in position 1 of
# dictionary value) with processed name text (element at 2nd position of list of dictionary value)
for DictionaryKey in DictionaryofIdsAndOldAndNewData.keys():
    modifiedfastalistfilestring = modifiedfastalistfilestring.replace((DictionaryofIdsAndOldAndNewData[DictionaryKey])[0],(DictionaryofIdsAndOldAndNewData[DictionaryKey])[1])
#for checking
#print modifiedfastalistfilestring[:4000]
 
 
 
#DEALING WITH DUPLIACTES
#Already removed key and value resulting from it in from dictionary tallies above
#Now remove any entries found to be duplicates from the modifiedfastalistfilestring
duplicateremovedmodifiedlistfilestring = "" #initialize emptey string to be new one
# I am just going to break it up at '>' at start of lines and put back together without including those flagged with the tag marking duplicates
ListOFEntriesInmodifiedfastalistfilestring = re.split ('\n>',modifiedfastalistfilestring) # use re.split because want it to match only beginning of lines to make more specific, just in case
del ListOFEntriesInmodifiedfastalistfilestring[0] #because the split is always going to create a list entry up ahead of first point the '>' occurs and then it will get an '>' added by next command.
ListOFEntriesInmodifiedfastalistfilestring = ["\n>"+ListOFEntriesInmodifiedfastalistfilestring[i] for i in range(len(ListOFEntriesInmodifiedfastalistfilestring)) ]  #see about 'LOCUS' above for why doing this way
#print ListOFEntriesInmodifiedfastalistfilestring
for eachentry in ListOFEntriesInmodifiedfastalistfilestring:
    if DuplicatesIndicatorTag not in eachentry:
        duplicateremovedmodifiedlistfilestring = duplicateremovedmodifiedlistfilestring + eachentry
#print duplicateremovedmodifiedlistfilestring[:4000]
sys.stderr.write("\nTotal entries that had different accessions but were duplicated organisms and were discarded: "+ str(NumberofSciNameDuplicatesNoted) +". ")
 
 
 
#Some more checks after merging names into main dictionary next to placeholder
#just a check that no values in the DictionaryofIdsAndOldAndNewData are still status 'NotLookedUp' when supposedly done processing list.
#SHOULD NOT TRIGGER.
if 'NotLookedUp' in (x[1] for x in DictionaryofIdsAndOldAndNewData.values()):    #from http://stackoverflow.com/questions/1156087/python-search-in-lists-of-lists
    sys.stderr.write("\nWARNING: Not all names were looked up and coding placeholders remain after lookup and merge.")
    SuggestReportBug()
if 'NotObtained' in (x[1] for x in DictionaryofIdsAndOldAndNewData.values()):    #from http://stackoverflow.com/questions/1156087/python-search-in-lists-of-lists
    sys.stderr.write("\nWARNING: Not all names were obtained during look up and coding placeholders remain.")
    SuggestReportBug()
 
 
#Following final replacement "PLACEHOLDER" (all caps) should not occur in the text
#SHOULD NOT TRIGGER.
if "PLACEHOLDER" in duplicateremovedmodifiedlistfilestring:
    sys.stderr.write("\nWARNING: Apparently not all entries have names added properly.")
    SuggestReportBug()
 
 
 
#give user some stats at end
"""
FinalScientificList = ScientificNameList + PDBScientificNameList
FinalOrganismNameList = OrganismNameList + PDBOrganismNameList
sys.stderr.write("\n\nScientific names related to entries: ")
for sciname in FinalScientificList:
    sys.stderr.write("     "+ sciname +"     ")
"""
 
 
 
 
 
 
#write results to a file
TheFileNameMainPart, fileExtension = os.path.splitext(FileName) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
if '.' in FileName:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
    OutputFileName=  TheFileNameMainPart+"TAXONED"+fileExtension
else:
    OutputFileName= FileName+"TAXONED"
modifiedfastafileoutput = open(rootpath + OutputFileName, "w")
modifiedfastafileoutput.write(duplicateremovedmodifiedlistfilestring)
sys.stderr.write("\nResults written to file named '"+ OutputFileName +"' ")
if rootpath != "":
    sys.stderr.write("at " + rootpath +" ")
sys.stderr.write(".")
modifiedfastafileoutput.close()
 
sys.stderr.write("\n")

"""

 
    Like to do:
    -  add Integrated Taxonomic Information System (http://www.itis.gov/) which
    # has an API described at http://www.itis.gov/web_service.html .
    I'd have to add use if that site if I wanted it 100% automated. 
    Presently PythonAnywhere doesn't have this on whitelist for free accounts.
 
 
 
 
"""
