#include <stdio.h>

#include "print_info.h"

void print_ELPH_logo(int mpi_rank, FILE* output)
{
    if (mpi_rank)
    {
        return;
    }
    fprintf(output,
            "       */*                                                        "
            "    *+        \n");
    fprintf(output,
            "      */                                                          "
            "     //       \n");
    fprintf(output,
            "    +/*   +*+            +          +++++++ ++   +++++   ++      "
            "+/***/  */+    \n");
    fprintf(output,
            "   */+    +*+     +++++ ***+ ++++*  ++      ++   ++   ++ ++++++  "
            "/*   +    /*   \n");
    fprintf(output,
            "  //      +*+    ** ++*  *+    +*   +++++++ ++   ++++++  ++  ++ "
            "+/+         */+ \n");
    fprintf(output,
            "  */+     +*+    **      *+   *+    ++      ++   ++      ++  ++  "
            "/*   *+    */* \n");
    fprintf(output,
            "   +/*     +++++  +++++  +++ +++++  +++++++ ++   ++      ++  ++   "
            "+***+    */   \n");
    fprintf(output,
            "     /*                                                           "
            "      +/*     \n");
    fprintf(output,
            "      */+            **//*+        **+     **        +**/*+       "
            "     */+      \n");
    fprintf(output,
            "       +/+         +//****//      /  /+   /  /      ////////+     "
            "    //        \n");
    fprintf(output,
            "                   +//+++*//+++++     /  /+  +*+*+++////////+     "
            "              \n");
    fprintf(output,
            "                    +/////*+          +**+           */////+      "
            "              \n");
    fprintf(output, "\n\n");
}
