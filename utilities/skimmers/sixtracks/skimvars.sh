grep -Eo 'fReader, "\w+' $1 > listofvars.txt #get variables
sed -i -e 's/fReader, \"//g' listofvars.txt
cp listofvars.txt listofvarsp4.txt
sed -i -e '/p4/d' listofvars.txt
sed -i -e '/p4/!d' listofvarsp4.txt

cp listofvars.txt first.txt
cp listofvars.txt mid.txt
cp listofvars.txt last.txt

cp listofvars.txt outvars.txt
cp listofvars.txt inpoint.txt

cp listofvarsp4.txt firstp4.txt
cp listofvarsp4.txt midp4.txt
cp listofvarsp4.txt lastp4.txt

cp listofvarsp4.txt outvarsp4.txt
cp listofvarsp4.txt inpointp4.txt

sed -i -e 's/^/out_/g' listofvars.txt #add out at beginning of each variables
sed -i -e 's/$/,/g' listofvars.txt #add , to separate vars
sed -i -e 's/^/out_/g' p4vars.txt #add out at beginning of each variables
sed -i -e 's/$/,/g' p4vars.txt #add , to separate vars

#Output variables
paste -d ' ' - - - - - < listofvars.txt > vars.txt #compact in one line each 5 lines
sed -i -e 's/^/Float_t /g' vars.txt #add Float_t ad the beginning of each line
sed -i -e 's/,$/;/g' vars.txt #replace , at the end of line with ;
sed -i -e '$ s/, $/;/g' vars.txt #replace , at the end of file with ;

paste -d ' ' - - - - - < p4vars.txt > varsp4.txt #compact in one line each 5 lines
sed -i -e 's/^/TLorentzVector /g' varsp4.txt #add Float_t ad the beginning of each line
sed -i -e 's/,$/;/g' varsp4.txt #replace , at the end of line with ;
sed -i -e '$ s/, $/;/g' varsp4.txt #replace , at the end of file with ;


#Output branches
sed -i -e 's/^/outTree->Branch("/g' first.txt && sed -i -e 's/$/", /g' first.txt
sed -i -e 's/^/\&out_/g' mid.txt && sed -i -e 's/$/, /g' mid.txt
sed -i -e 's/^/\"/g' last.txt && sed -i -e 's/$/\/F\");/g' last.txt

#dimuon_tree->Branch("lowMuon_p4",  "TLorentzVector", &lowMuon_p4);
sed -i -e 's/^/outTree->Branch("/g' firstp4.txt && sed -i -e 's/$/", "TLorentzVector", /g' firstp4.txt
sed -i -e 's/^/\&/g' lastp4.txt && sed -i -e 's/$/\");/g' lastp4.txt

paste firstp4.txt lastp4.txt > branchesp4.txt

#Variables
sed -i -e 's/^/out_/g' outvars.txt  && sed -i -e 's/$/ = /g' outvars.txt
sed -i -e 's/^/(Float_t)(*/g' inpoint.txt  && sed -i -e 's/$/);/g' inpoint.txt

paste outvars.txt inpoint.txt > assign.txt

sed -i -e 's/^/out_/g' outvarsp4.txt  && sed -i -e 's/$/ = /g' outvarsp4.txt
sed -i -e 's/^/(*/g' inpointp4.txt  && sed -i -e 's/$/);/g' inpointp4.txt

paste outvarsp4.txt inpointp4.txt > assignp4.txt

#Cleaning
rm first*.txt
rm mid*.txt
rm last*.txt
rm listofvars*.txt
rm outvars*.txt
rm inpoint*.txt
