#pragma once

// LC_NUMERIC switched to a locale with a decimal comma, for the tests of what must not depend on it. The
// process locale can change under mdlib: a GTK file dialog switches it to the user's, for one.

#include <locale.h>
#include <stdbool.h>
#include <stdio.h>

typedef struct comma_locale_t {
    char prev[256];
} comma_locale_t;

// False, with nothing changed, when no decimal comma locale is installed
static inline bool comma_locale_enter(comma_locale_t* loc) {
    const char* prev = setlocale(LC_NUMERIC, NULL);
    snprintf(loc->prev, sizeof(loc->prev), "%s", prev ? prev : "C");

    static const char* names[] = { "sv_SE.UTF-8", "sv_SE.utf8", "de_DE.UTF-8", "de_DE.utf8", "sv-SE", "de-DE", "Swedish_Sweden.1252" };
    for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); ++i) {
        if (setlocale(LC_NUMERIC, names[i]) && localeconv()->decimal_point[0] == ',') {
            return true;
        }
    }
    setlocale(LC_NUMERIC, loc->prev);
    printf("  no decimal comma locale installed, nothing to check\n");
    return false;
}

static inline void comma_locale_leave(const comma_locale_t* loc) {
    setlocale(LC_NUMERIC, loc->prev);
}
