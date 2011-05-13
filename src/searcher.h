/* 
 * File:   searcher.h
 * Author: david
 *
 * Created on 13 mei 2011, 10:04
 */

#ifndef SEARCHER_H
#define	SEARCHER_H

bool CheckMismatches(const std::string & TheInput,
        const std::string & CurrentReadSeq,
        const UniquePoint & UP);

void CheckLeft_Far(SPLIT_READ & OneRead,
        const std::string & TheInput,
        const std::string & CurrentReadSeq,
        const std::vector <unsigned int> Left_PD[],
        const short & BP_Left_Start,
        const short & BP_Left_End,
        const short & CurrentLength,
        std::vector <UniquePoint> & LeftUP);

void CheckRight_Far(SPLIT_READ & OneRead,
        const std::string & TheInput,
        const std::string & CurrentReadSeq,
        const std::vector <unsigned int> Right_PD[],
        const short & BP_Right_Start,
        const short & BP_Right_End,
        const short & CurrentPos,
        std::vector <UniquePoint> & RightUP);

void CheckLeft_Close(const SPLIT_READ & OneRead,
        const std::string & TheInput,
        const std::string & CurrentReadSeq,
        const std::vector <unsigned int> Left_PD[],
        const short & BP_Left_Start,
        const short & BP_Left_End,
        const short & CurrentLength,
        std::vector <UniquePoint> & LeftUP);

void CheckRight_Close(const SPLIT_READ & OneRead,
        const std::string & TheInput,
        const std::string & CurrentReadSeq,
        const std::vector <unsigned int> Right_PD[],
        const short & BP_Right_Start,
        const short & BP_Right_End,
        const short & CurrentPos,
        std::vector <UniquePoint> & RightUP);

#endif	/* SEARCHER_H */

