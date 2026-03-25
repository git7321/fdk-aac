/******************* Library for basic calculation routines ********************

   Author(s):   Manuel Jander

   Description: FDK tools versioning support

*******************************************************************************/

#include "FDK_core.h"

#define FDK_TOOLS_LIB_VL0 3
#define FDK_TOOLS_LIB_VL1 1
#define FDK_TOOLS_LIB_VL2 0
#define FDK_TOOLS_LIB_TITLE "FDK Tools"
#define FDK_TOOLS_LIB_BUILD_DATE ""
#define FDK_TOOLS_LIB_BUILD_TIME ""

int FDK_toolsGetLibInfo(LIB_INFO *info) {
  UINT v;
  int i;

  if (info == NULL) {
    return -1;
  }

  i = FDKlibInfo_lookup(info, FDK_TOOLS);
  if (i < 0) return -1;

  info += i;

  v = LIB_VERSION(FDK_TOOLS_LIB_VL0, FDK_TOOLS_LIB_VL1, FDK_TOOLS_LIB_VL2);

  FDKsprintf(info->versionStr, "%d.%d.%d", ((v >> 24) & 0xff),
             ((v >> 16) & 0xff), ((v >> 8) & 0xff));

  info->module_id = FDK_TOOLS;
  info->version = v;
  info->build_date = FDK_TOOLS_LIB_BUILD_DATE;
  info->build_time = FDK_TOOLS_LIB_BUILD_TIME;
  info->title = FDK_TOOLS_LIB_TITLE;
  info->flags = 1;

  return 0;
}
#if defined(_WIN32) && defined(FDK_EXPORT)
#include <windows.h>

BOOL APIENTRY DllMain(HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved)
{
	switch (ul_reason_for_call)
	{
	    case DLL_PROCESS_ATTACH:
	    case DLL_THREAD_ATTACH:
	    case DLL_THREAD_DETACH:
	    case DLL_PROCESS_DETACH:
		    break;
	}
	return TRUE;
}

#endif
