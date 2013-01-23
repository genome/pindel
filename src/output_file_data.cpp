#include <iostream>

#include "output_file_data.h"

OutputFileData::OutputFileData()
{
   m_templateCounter = 0;
   m_nonTemplateCounter = 0;
}

int OutputFileData::getSvIndex() const
{
   return m_templateCounter + m_nonTemplateCounter;
}

int OutputFileData::getTemplateSvCounter() const
{
   return m_templateCounter;
}

int OutputFileData::getNonTemplateSvCounter() const
{
   return m_nonTemplateCounter;
}

void OutputFileData::increaseTemplateSvCounter()
{
   m_templateCounter++;
}

void OutputFileData::increaseNonTemplateSvCounter()
{
   m_nonTemplateCounter++;
}
