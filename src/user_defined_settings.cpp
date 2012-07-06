/* UserDefinedSettings' manages the parameters that have been set by the users or have predefined default values. */

#include "user_defined_settings.h"

UserDefinedSettings* UserDefinedSettings::m_instance = NULL;

UserDefinedSettings::UserDefinedSettings()
{
}

UserDefinedSettings* UserDefinedSettings::Instance() {
	if (m_instance == NULL ) {
		m_instance = new UserDefinedSettings;
	}
	return m_instance;
}
