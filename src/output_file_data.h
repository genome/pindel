/* For every output file, save data on the counters so that counters are unique and increased correctly.

	Created by Eric-Wubbo Lameijer, Leiden University Medical Center, July 21st 2011, e.m.w.lameijer@lumc.nl */

class OutputFileData {
public:
	OutputFileData();
   int getSvIndex() const;
	int getTemplateSvCounter() const;
	int getNonTemplateSvCounter() const;
	void increaseTemplateSvCounter();
	void increaseNonTemplateSvCounter();

private:
   int m_templateCounter;
	int m_nonTemplateCounter;
};




