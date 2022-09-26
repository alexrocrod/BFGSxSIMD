// #include "stdafx.h"
#define AADC_INTERNAL_DEBUG


//#include "aadc/aadc_debug.h"
#include <aadc/aadc_matrix.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <immintrin.h>
#include <unordered_map>

#include <memory>
#include <cstdlib>

#include <aadc/idouble.h>
#include <aadc/ibool.h>
#include <aadc/iint.h>

#include <random>

#include <aadc/aadc.h>
// #include <aadc/aadc_private.h>
#include <ql/quotes/simplequote.hpp>
#include <ql/quantlib.hpp>

#include <ql/time/daycounter.hpp>
#include <ql/time/daycounters/all.hpp>
#include <ql/time/calendars/jointcalendar.hpp>
#include <ql/time/calendars/all.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>

#include "CurveCalibrate.h"

namespace QuantLib {
    namespace AAD {
        typedef QuantLib::Calendar Calendar;
        typedef QuantLib::UnitedKingdom UnitedKingdom;
        typedef QuantLib::UnitedStates UnitedStates;
    }
}

int main(int argc, char* argv[]) {
using namespace QuantLib;

    Calendar calendar = JointCalendar(UnitedKingdom(UnitedKingdom::Exchange),UnitedStates(UnitedStates::Settlement), JoinHolidays);
	Date settlementDate(18, February, 2014);
	settlementDate = calendar.adjust(settlementDate);
	Integer fixingDays = 2;
	Date todaysDate = calendar.advance(settlementDate, -fixingDays, Days);
	Settings::instance().evaluationDate() = todaysDate;
	DayCounter depositDayCounter = Actual360();

	typedef __m256d mmType;
	DualOISCurveCalibration calibration(settlementDate);
	calibration.recordPortfolioRevalue();
	calibration.calibrateCurveNewton();
	for (int i=0; i<7; i++) {
		std::cout << "Exp " << i <<"\n";
		calibration.PricingAfterCalibration(i);
	}

	return 0; // runMain<__m256d>();
}

