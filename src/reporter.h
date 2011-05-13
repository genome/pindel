/* 
 * File:   reporter.h
 * Author: david
 *
 * Created on 13 mei 2011, 11:27
 */

#ifndef REPORTER_H
#define	REPORTER_H

void SortOutputD(const unsigned & NumBoxes, const std::string & CurrentChr, std::vector <SPLIT_READ> & AllReads, std::vector <unsigned> Deletions[], std::ofstream & DeletionOutf);
void SortOutputSI(const unsigned & NumBoxes, const std::string & CurrentChr, std::vector <SPLIT_READ> & AllReads, std::vector <unsigned> SIs[], std::ofstream & SIsOutf);
void SortOutputTD(const unsigned & NumBoxes, const std::string & CurrentChr, std::vector <SPLIT_READ> & AllReads, std::vector <unsigned> TDs[], std::ofstream & TDOutf);
void SortOutputTD_NT(const unsigned & NumBoxes, const std::string & CurrentChr, std::vector <SPLIT_READ> & AllReads, std::vector <unsigned> TDs[], std::ofstream & TDOutf);
void SortOutputInv(const unsigned & NumBoxes, const std::string & CurrentChr, std::vector <SPLIT_READ> & AllReads, std::vector <unsigned> Inv[], std::ofstream & InvOutf);
void SortOutputInv_NT(const unsigned & NumBoxes, const std::string & CurrentChr, std::vector <SPLIT_READ> & AllReads, std::vector <unsigned> Inv[], std::ofstream & InvOutf);
void SortOutputDI(const unsigned & NumBoxes, const std::string & CurrentChr, std::vector <SPLIT_READ> & AllReads, std::vector <unsigned> DI[], std::ofstream & DIOutf);
void SortOutputLI(const std::string & CurrentChr, std::vector <SPLIT_READ> & AllReads, std::ofstream & Outf_LI);
void SortOutputRest(const std::string & CurrentChr, std::vector <SPLIT_READ> & AllReads, std::vector <SPLIT_READ> & BP_Reads, std::ofstream & Outf_Rest);


#endif	/* REPORTER_H */

