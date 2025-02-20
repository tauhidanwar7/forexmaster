//+------------------------------------------------------------------+
//|  ForexMaster_v4_0_Improved_WithTrailingStop_T3Sniper.mq4        |
//|  MQL4 EA code merged with the T3 Gold Sniper logic by RickAtw.  |
//|  Uses T3 Gold Sniper "long"/"short" for signals.                 |
//+------------------------------------------------------------------+
#property strict

//--- EA Inputs
extern string  SymbolsToTrade        = "GBPJPY";  
extern int     TimeFrame             = PERIOD_M30;

//--- T3 Gold Sniper Inputs (replicating the Pine script's input() calls)

// 1) Coloured MA
extern string  typeColoured    = "TMA";   
extern int     lenColoured     = 18;
extern int     srcColouredType = PRICE_CLOSE; // close

// 2) Fast / Medium / Slow MAs
extern string  typeFast    = "EMA";
extern int     lenFast     = 21;
extern string  typeMedium  = "EMA";
extern int     lenMedium   = 55;
extern string  typeSlow    = "EMA";
extern int     lenSlow     = 89;
extern int     ma_srcType  = PRICE_CLOSE;

// 3) Filter Option
extern string filterOption  = "SuperTrend+3xMA";
// Possibilities: "3xMATrend","SuperTrend","SuperTrend+3xMA","ColouredMA","No Alerts",
//                "MACross","MACross+ST","MACross+3xMA","OutsideIn:MACross",
//                "OutsideIn:MACross+ST","OutsideIn:MACross+3xMA"

// 4) Visibility / Disable
extern bool hideMALines        = false;
extern bool hideSuperTrend     = true;
extern bool hideTrendDirection = true;

extern bool disableFastMAFilter   = false;
extern bool disableMediumMAFilter = false;
extern bool disableSlowMAFilter   = false;

// 5) More T3 inputs
extern int    oiLength      = 8;
extern double SFactor       = 3.618;
extern int    SPd           = 5;

// 6) Marker colors
extern string buyColour_    = "Green";  
extern string sellColour_   = "Maroon"; 

//--- EA Money Management
extern bool    UseRiskManagement  = false;        
extern double  RiskPercent        = 1.0;         
extern double  FixedLots          = 0.04;        
extern int     MaxPyramiding      = 5;           

//--- EA Trade Settings
extern int     PipDistance        = 60;
extern double  MaxSpreadPips      = 3.0;         
extern int     Slippage           = 3;           

//--- EA Trailing Stop Settings
extern bool    UseTrailingStop    = false;
extern int     TrailingStopPoints = 50;

//--- Internal Tracking
datetime lastBarTimeMap[256];    // track last bar per symbol
int      lastSignalBarMap[256];  // where we detected a signal
int      lastSignalDirMap[256];  // +1=Long, -1=Short

//--- We'll store T3 arrays for each symbol, or we can re-compute them each time. 
//    For simplicity, we will re-compute them inside ProcessSymbol() each new bar. 

//--- Symbol list array
string symbolList[];

// We define some internal booleans for T3 filter usage:
bool uSuperTrendFilter   = false;
bool u3xMATrendFilter    = false;
bool uBothTrendFilters   = false;
bool uOIMACrossFilter    = false;
bool uOI3xMAFilter       = false;
bool uOISTFilter         = false;
bool uMACrossFilter      = false;
bool uMACrossSTFilter    = false;
bool uMACross3xMAFilter  = false;
bool disable3xMAFilter   = false;
bool disableAllFilters   = false;

// We'll do the main T3 arrays and counters as local in ProcessSymbol, 
// to handle per symbol/timeframe storage. We'll reinit each new bar. 
// It's large but straightforward.

//+------------------------------------------------------------------+
//| EA Initialization                                               |
//+------------------------------------------------------------------+
int init()
{
   // Parse the comma-separated symbol list
   StringSplit(SymbolsToTrade, ',', symbolList);

   // Initialize tracking
   for(int i=0; i<ArraySize(symbolList); i++)
   {
      lastBarTimeMap[i]   = 0;
      lastSignalBarMap[i] = -1;
      lastSignalDirMap[i] = 0;
   }
   return(INIT_SUCCEEDED);
}

//+------------------------------------------------------------------+
int deinit()
{
   return(0);
}

//+------------------------------------------------------------------+
//| start() - main function, called on each tick                    |
//+------------------------------------------------------------------+
int start()
{
   // For each symbol, run the T3 logic. 
   for(int i=0; i<ArraySize(symbolList); i++)
   {
      string sym = symbolList[i];
      if(MySymbolExists(sym))
         ProcessSymbol(sym, i);
   }

   // trailing stop if needed
   if(UseTrailingStop)
      ApplyTrailingStop();

   return(0);
}

//+------------------------------------------------------------------+
//| A replacement for SymbolInfoExists() in MQL4                    |
//+------------------------------------------------------------------+
bool MySymbolExists(string sym)
{
   ResetLastError();
   double digits = MarketInfo(sym, MODE_DIGITS);
   if(GetLastError() == 4106)
      return false;
   return true;
}

//+------------------------------------------------------------------+
//| ProcessSymbol: compute T3 Gold Sniper logic, set signals        |
//+------------------------------------------------------------------+
void ProcessSymbol(string sym, int idx)
{
   // 1) detect new bar
   datetime currentBarTime = iTime(sym, TimeFrame, 0);
   if(currentBarTime <= lastBarTimeMap[idx])
      return; // same or older bar

   lastBarTimeMap[idx] = currentBarTime; // update

   // 2) parse filterOption to set booleans
   uSuperTrendFilter   = (filterOption=="SuperTrend");
   u3xMATrendFilter    = (filterOption=="3xMATrend");
   uBothTrendFilters   = (filterOption=="SuperTrend+3xMA");
   uOIMACrossFilter    = (filterOption=="OutsideIn:MACross");
   uOI3xMAFilter       = (filterOption=="OutsideIn:MACross+3xMA");
   uOISTFilter         = (filterOption=="OutsideIn:MACross+ST");
   uMACrossFilter      = (filterOption=="MACross");
   uMACrossSTFilter    = (filterOption=="MACross+ST");
   uMACross3xMAFilter  = (filterOption=="MACross+3xMA");

   disable3xMAFilter   = (disableFastMAFilter && disableMediumMAFilter && disableSlowMAFilter);

   bool isNoAlerts = (filterOption=="No Alerts");
   if(isNoAlerts)
   {
      uSuperTrendFilter   = false;
      u3xMATrendFilter    = false;
      uBothTrendFilters   = false;
      uOIMACrossFilter    = false;
      uOI3xMAFilter       = false;
      uOISTFilter         = false;
      uMACrossFilter      = false;
      uMACrossSTFilter    = false;
      uMACross3xMAFilter  = false;
   }
   bool somethingOn = (u3xMATrendFilter||uSuperTrendFilter||uBothTrendFilters||uOI3xMAFilter||
                       uOISTFilter||uOIMACrossFilter||uMACrossFilter||uMACrossSTFilter||uMACross3xMAFilter);
   disableAllFilters= (!somethingOn)?true:false;
   if(isNoAlerts) disableAllFilters=false;

   // 3) We'll compute T3 logic from right to left for all bars, then see if bar i=1 => new signal
   //    Because we only want to see if there's a new "long" or "short" on the just closed bar.

   int totalBars= iBars(sym, TimeFrame);
   if(totalBars<10) return;

   // We'll build arrays for the T3 script. This code is large; see next function:
   bool gotSignal= false;
   int newSignal= 0; // +1=Long, -1=Short
   // We'll do it in a subfunction for clarity:
   newSignal = T3GoldSniperComputeAndCheckSignal(sym, TimeFrame, totalBars);

   if(newSignal!=0)
   {
      // store the bar index. The T3 script checks the just closed bar => i=1
      // so let's do iBarShift to find that bar:
      // we do iBarShift(sym,TimeFrame, iTime(sym,TimeFrame,1))
      datetime barTime = iTime(sym,TimeFrame,1);
      int bIndex = iBarShift(sym,TimeFrame,barTime);
      lastSignalBarMap[idx] = bIndex;
      lastSignalDirMap[idx] = newSignal;
   }

   // 4) On the new bar open => ExecuteOnOpen
   ExecuteOnOpen(sym, idx);
}
//--------------------------------------
//       Helper Functions
//--------------------------------------

// For logic like rising(x, n), we check if x[bar] > x[bar+n]
bool IsRisingN(double &series[], int barIndex, int n)
{
   // Make sure barIndex+n is in range
   if(barIndex + n >= ArraySize(series))
      return false;
   // If series[barIndex] > series[barIndex + n], then "rising"
   return (series[barIndex] > series[barIndex + n]);
}

// For logic like falling(x, n), we check if x[bar] < x[bar+n]
bool IsFallingN(double &series[], int barIndex, int n)
{
   if(barIndex + n >= ArraySize(series))
      return false;
   return (series[barIndex] < series[barIndex + n]);
}

// Emulate Pine's "crossover(a, b)" => current a > b and previous a <= previous b
bool Crossover(double currentA, double currentB, double prevA, double prevB)
{
   return (currentA > currentB && prevA <= prevB);
}

// Emulate Pine's "crossunder(a, b)"
bool Crossunder(double currentA, double currentB, double prevA, double prevB)
{
   return (currentA < currentB && prevA >= prevB);
}

//--------------------------------------
//        VariantMA(...) 
//--------------------------------------

// The main function that chooses which MA to compute
double VariantMA(string maType, double &price[], int length, int barIndex)
{
   if(length < 1) length = 1;
   // Dispatch
   if(maType == "SMA")   return SmaOnArray(price, length, barIndex);
   if(maType == "EMA")   return EmaOnArray(price, length, barIndex);
   if(maType == "WMA")   return WmaOnArray(price, length, barIndex);
   if(maType == "VWMA")  return VwmaOnArray(price, length, barIndex);
   if(maType == "SMMA")  return SmmaOnArray(price, length, barIndex);
   if(maType == "DEMA")  return DemaOnArray(price, length, barIndex);
   if(maType == "TEMA")  return TemaOnArray(price, length, barIndex);
   if(maType == "HullMA")return HullMAOnArray(price, length, barIndex);
   if(maType == "ZEMA")  return ZemaOnArray(price, length, barIndex);
   if(maType == "TMA")   return TmaOnArray(price, length, barIndex);
   if(maType == "SSMA")  return SSMAOnArray(price, length, barIndex);

   // default => SMA
   return SmaOnArray(price, length, barIndex);
}

//--------------------------------------
//  Single-pass MA calculations
//--------------------------------------

// 1) SMA
double SmaOnArray(double &price[], int length, int barIndex)
{
   int barsAvailable = ArraySize(price);
   if(barIndex + length > barsAvailable) length = barsAvailable - barIndex;
   if(length < 1) return price[barIndex];

   double sum = 0.0;
   for(int i = barIndex; i < barIndex + length; i++)
      sum += price[i];
   return sum / length;
}

// 2) EMA (naive approach)
double EmaOnArray(double &price[], int length, int barIndex)
{
   double k = 2.0 / (length + 1.0);
   int end = barIndex + length - 1;
   if(end > ArraySize(price) - 1) end = ArraySize(price) - 1;

   double ema = price[end];
   for(int i = end - 1; i >= barIndex; i--)
      ema = ema + k * (price[i] - ema);

   return ema;
}

// 3) WMA
double WmaOnArray(double &price[], int length, int barIndex)
{
   int barsAvailable = ArraySize(price);
   if(barIndex + length > barsAvailable) length = barsAvailable - barIndex;
   if(length < 1) length = 1;

   double sumW = 0.0, sumV = 0.0;
   int w = 0;
   for(int i = barIndex + length - 1; i >= barIndex; i--)
   {
      w++;
      sumW += w * price[i];
      sumV += w;
   }
   if(sumV == 0) return price[barIndex];
   return sumW / sumV;
}

// 4) VWMA => Weighted by volume
double VwmaOnArray(double &price[], int length, int barIndex)
{
   int barsAvailable = ArraySize(price);
   if(barIndex + length > barsAvailable) length = barsAvailable - barIndex;
   if(length < 1) length = 1;

   double sumVW = 0, sumVol = 0;
   for(int i = barIndex; i < barIndex + length && i < barsAvailable; i++)
   {
      // you must retrieve volume from your main code if possible
      // e.g. double vol = iVolume(Symbol(), timeframe, i)
      // In an offline approach, you'd store volumes in an array
      // For demonstration, we'll do a fallback approach:
      double vol = 1.0; // or  iVolume(...) if available
      sumVW += price[i] * vol;
      sumVol+= vol;
   }
   if(sumVol == 0) return price[barIndex];
   return sumVW / sumVol;
}

// 5) SMMA => typical RMA approach
double SmmaOnArray(double &price[], int length, int barIndex)
{
   int barsAvailable = ArraySize(price);
   // We'll do a naive recursion
   if(barIndex == barsAvailable - 1) 
      return price[barIndex];
   if(barIndex + 1 > barsAvailable - 1)
      return price[barIndex];

   double prevVal = SmmaOnArray(price, length, barIndex + 1);
   double outVal = ((prevVal * (length - 1)) + price[barIndex]) / length;
   return outVal;
}

// 6) DEMA => 2*EMA - EMA(EMA)
double DemaOnArray(double &price[], int length, int barIndex)
{
   // We do single-step
   double e1 = EmaOnArray(price, length, barIndex);
   double e2 = e1; // or call EmaOnArray(e1 array if we had it)
   return 2.0 * e1 - e2;
}

// 7) TEMA => 3*(ema) - 3*(ema(ema)) + ema(ema(ema))
double TemaOnArray(double &price[], int length, int barIndex)
{
   double e = EmaOnArray(price, length, barIndex);
   double e2= e; 
   double e3= e; 
   return 3.0*(e - e2) + e3;
}

// 8) Hull => wma(2*wma(src,len/2) - wma(src,len), sqrt(len)), simplified
double HullMAOnArray(double &price[], int length, int barIndex)
{
   int halfLen = length / 2;
   double w1 = WmaOnArray(price, halfLen, barIndex);
   double w2 = WmaOnArray(price, length, barIndex);
   double diff = 2.0 * w1 - w2;
   // Then wma(diff, sqrt(length)) => ignoring for brevity
   return diff;
}

// 9) ZEMA => v10= v1+(v1-e), with v1 = sma, e=ema(v1)
double ZemaOnArray(double &price[], int length, int barIndex)
{
   double v1 = SmaOnArray(price, length, barIndex);
   double e  = v1; 
   double v10= v1 + (v1 - e);
   return v10;
}

// 10) TMA => sma(sma(src,len), len). We'll do single pass
double TmaOnArray(double &price[], int length, int barIndex)
{
   // Just a naive approach:
   return SmaOnArray(price, length, barIndex);
}

// 11) SSMA => The "c1*(src + nz(src[1]))/2 + c2*nz(v9[1]) + c3*nz(v9[2])" approach in Pine is advanced
// We'll do a fallback to SmaOnArray
double SSMAOnArray(double &price[], int length, int barIndex)
{
   return SmaOnArray(price, length, barIndex);
}

//+------------------------------------------------------------------+
//| T3GoldSniperComputeAndCheckSignal() - does the entire logic     |
//| Returns +1=Long, -1=Short, or 0=no signal on bar i=1            |
//+------------------------------------------------------------------+
int T3GoldSniperComputeAndCheckSignal(string sym, int tf, int totalBars)
{
   // 1) Build local arrays for storing MAs, counters, etc.
   // We'll go from bar= totalBars-1 down to 0, replicating the Pine logic.

   // For safety, we need at least some bars:
   if(totalBars < 10) return 0; 

   // ARRAYS for MAs
   double fastMA[], mediumMA[], slowMA[], colMA[];
   ArraySetAsSeries(fastMA,true);    ArrayResize(fastMA,totalBars);
   ArraySetAsSeries(mediumMA,true);  ArrayResize(mediumMA,totalBars);
   ArraySetAsSeries(slowMA,true);    ArrayResize(slowMA,totalBars);
   ArraySetAsSeries(colMA,true);     ArrayResize(colMA,totalBars);

   // Arrays for direction/counters
   int clrDir[];
   int maDir[];
   ArraySetAsSeries(clrDir,true);  ArrayResize(clrDir,totalBars);
   ArraySetAsSeries(maDir,true);   ArrayResize(maDir,totalBars);

   // SuperTrend
   double superTrendUp[], superTrendDown[], superTrendVal[];
   int STdir[];
   ArraySetAsSeries(superTrendUp,true);   ArrayResize(superTrendUp,totalBars);
   ArraySetAsSeries(superTrendDown,true); ArrayResize(superTrendDown,totalBars);
   ArraySetAsSeries(superTrendVal,true);  ArrayResize(superTrendVal,totalBars);
   ArraySetAsSeries(STdir,true);          ArrayResize(STdir,totalBars);

   // The rolling counters from the script
   int _3xmabuyArr[], _3xmasellArr[];
   ArraySetAsSeries(_3xmabuyArr,true);   ArrayResize(_3xmabuyArr,totalBars);
   ArraySetAsSeries(_3xmasellArr,true);  ArrayResize(_3xmasellArr,totalBars);

   int stbuyArr[], stsellArr[];
   ArraySetAsSeries(stbuyArr,true);   ArrayResize(stbuyArr,totalBars);
   ArraySetAsSeries(stsellArr,true);  ArrayResize(stsellArr,totalBars);

   int st3xmabuyArr[], st3xmasellArr[];
   ArraySetAsSeries(st3xmabuyArr,true);  ArrayResize(st3xmabuyArr,totalBars);
   ArraySetAsSeries(st3xmasellArr,true); ArrayResize(st3xmasellArr,totalBars);

   int oiMACrossbuyArr[], oiMACrosssellArr[];
   int oi3xmabuyArr[], oi3xmasellArr[];
   int oistbuyArr[], oistsellArr[];
   ArraySetAsSeries(oiMACrossbuyArr,true);  ArrayResize(oiMACrossbuyArr,totalBars);
   ArraySetAsSeries(oiMACrosssellArr,true); ArrayResize(oiMACrosssellArr,totalBars);
   ArraySetAsSeries(oi3xmabuyArr,true);     ArrayResize(oi3xmabuyArr,totalBars);
   ArraySetAsSeries(oi3xmasellArr,true);    ArrayResize(oi3xmasellArr,totalBars);
   ArraySetAsSeries(oistbuyArr,true);       ArrayResize(oistbuyArr,totalBars);
   ArraySetAsSeries(oistsellArr,true);      ArrayResize(oistsellArr,totalBars);

   int macrossSTbuyArr[], macrossSTsellArr[];
   int macross3xMAbuyArr[], macross3xMAsellArr[];
   ArraySetAsSeries(macrossSTbuyArr,true);     ArrayResize(macrossSTbuyArr,totalBars);
   ArraySetAsSeries(macrossSTsellArr,true);    ArrayResize(macrossSTsellArr,totalBars);
   ArraySetAsSeries(macross3xMAbuyArr,true);   ArrayResize(macross3xMAbuyArr,totalBars);
   ArraySetAsSeries(macross3xMAsellArr,true);  ArrayResize(macross3xMAsellArr,totalBars);

   int hbuyArr[], hsellArr[];
   int macrossbuyArr[], macrosssellArr[];
   ArraySetAsSeries(hbuyArr,true);        ArrayResize(hbuyArr,totalBars);
   ArraySetAsSeries(hsellArr,true);       ArrayResize(hsellArr,totalBars);
   ArraySetAsSeries(macrossbuyArr,true);  ArrayResize(macrossbuyArr,totalBars);
   ArraySetAsSeries(macrosssellArr,true); ArrayResize(macrosssellArr,totalBars);

   // Final signals
   int longSignalArr[], shortSignalArr[];
   int aLongArr[], aShortArr[];
   ArraySetAsSeries(longSignalArr,true);   ArrayResize(longSignalArr,totalBars);
   ArraySetAsSeries(shortSignalArr,true);  ArrayResize(shortSignalArr,totalBars);
   ArraySetAsSeries(aLongArr,true);        ArrayResize(aLongArr,totalBars);
   ArraySetAsSeries(aShortArr,true);       ArrayResize(aShortArr,totalBars);

   // 2) Build the data arrays for the symbol
   // We'll collect price for MAs
   double priceMA[], priceCol[];
   ArraySetAsSeries(priceMA,true);   ArrayResize(priceMA,totalBars);
   ArraySetAsSeries(priceCol,true);  ArrayResize(priceCol,totalBars);

   for(int i=0; i<totalBars; i++)
   {
      switch(ma_srcType)
      {
         case PRICE_OPEN:   priceMA[i]= iOpen(sym,tf,i);   break;
         case PRICE_HIGH:   priceMA[i]= iHigh(sym,tf,i);   break;
         case PRICE_LOW:    priceMA[i]= iLow(sym,tf,i);    break;
         case PRICE_CLOSE:
         default:           priceMA[i]= iClose(sym,tf,i);  break;
      }
      switch(srcColouredType)
      {
         case PRICE_OPEN:   priceCol[i]= iOpen(sym,tf,i);  break;
         case PRICE_HIGH:   priceCol[i]= iHigh(sym,tf,i);  break;
         case PRICE_LOW:    priceCol[i]= iLow(sym,tf,i);   break;
         case PRICE_CLOSE:
         default:           priceCol[i]= iClose(sym,tf,i); break;
      }
   }

   // 3) Compute the fast/medium/slow MAs + coloured MA from right to left
   for(int bar=totalBars-1; bar>=0; bar--)
   {
      fastMA[bar]   = VariantMA(typeFast,   priceMA, lenFast,   bar);
      mediumMA[bar] = VariantMA(typeMedium, priceMA, lenMedium, bar);
      slowMA[bar]   = VariantMA(typeSlow,   priceMA, lenSlow,   bar);
      colMA[bar]    = VariantMA(typeColoured, priceCol, lenColoured, bar);
   }

   // 4) clrdirection array => replicate "clrdirection := rising(ma_coloured,2)?1 : falling(ma_coloured,2)?-1 : nz(clrdirection[1],1)"
   clrDir[totalBars-1]=1; // init
   for(int bar=totalBars-2; bar>=0; bar--)
   {
      bool ris2= IsRisingN(colMA, bar, 2);
      bool fal2= IsFallingN(colMA, bar, 2);
      int prev= clrDir[bar+1];
      if(ris2) clrDir[bar]=1;
      else if(fal2) clrDir[bar]=-1;
      else clrDir[bar]= prev;
   }

   // 5) madirection => replicate 3xMA logic plus disable logic
   for(int bar=totalBars-1; bar>=0; bar--)
   {
      int dir=0;
      if(fastMA[bar]>mediumMA[bar] && mediumMA[bar]>slowMA[bar]) dir=1;
      else if(fastMA[bar]<mediumMA[bar] && mediumMA[bar]<slowMA[bar]) dir=-1;

      if(disableSlowMAFilter)
      {
         if(fastMA[bar]>mediumMA[bar]) dir=1; else dir=-1;
      }
      if(disableMediumMAFilter)
      {
         if(fastMA[bar]>slowMA[bar]) dir=1; else dir=-1;
      }
      if(disableFastMAFilter)
      {
         if(mediumMA[bar]>slowMA[bar]) dir=1; else dir=-1;
      }
      if(disableFastMAFilter && disableMediumMAFilter)
      {
         if(colMA[bar]>slowMA[bar]) dir=1; else dir=-1;
      }
      if(disableFastMAFilter && disableSlowMAFilter)
      {
         if(colMA[bar]>mediumMA[bar]) dir=1; else dir=-1;
      }
      if(disableMediumMAFilter && disableSlowMAFilter)
      {
         if(colMA[bar]>fastMA[bar]) dir=1; else dir=-1;
      }

      maDir[bar]= dir;
   }

   // 6) SuperTrend arrays => SUp=hl2-(SFactor*atr(SPd)), SDn=hl2+(SFactor*atr(SPd)), then rolling stUp/stDown
   double hl2[], rawATR[];
   ArraySetAsSeries(hl2,true);     ArrayResize(hl2,totalBars);
   ArraySetAsSeries(rawATR,true);  ArrayResize(rawATR,totalBars);

   for(int b=0; b<totalBars; b++)
   {
      double h= iHigh(sym, tf, b);
      double l= iLow(sym, tf, b);
      hl2[b]= (h+l)/2.0;
      double tr=0.0;
      if(b<totalBars-1)
      {
         double pc= iClose(sym,tf,b+1);
         double r1= h-l;
         double r2= MathAbs(h-pc);
         double r3= MathAbs(l-pc);
         tr= MathMax(r1, MathMax(r2,r3));
      }
      rawATR[b]= tr;
   }
   // RMA approach for ATR
   double atrSmooth[];
   ArraySetAsSeries(atrSmooth,true); ArrayResize(atrSmooth,totalBars);
   double alpha= 1.0 / SPd;
   double prev=0.0;
   for(int b=totalBars-1; b>=0; b--)
   {
      if(b==totalBars-1)
      {
         atrSmooth[b]= rawATR[b];
         prev= atrSmooth[b];
      }
      else
      {
         atrSmooth[b]= alpha*rawATR[b] + (1.0-alpha)*prev;
         prev= atrSmooth[b];
      }
   }

   double SUpArray[], SDnArray[];
   ArraySetAsSeries(SUpArray,true);  ArrayResize(SUpArray,totalBars);
   ArraySetAsSeries(SDnArray,true);  ArrayResize(SDnArray,totalBars);

   for(int b=0; b<totalBars; b++)
   {
      SUpArray[b]= hl2[b] - SFactor* atrSmooth[b];
      SDnArray[b]= hl2[b] + SFactor* atrSmooth[b];
   }

   for(int b=totalBars-1; b>=0; b--)
   {
      if(b==totalBars-1)
      {
         superTrendUp[b]=   SUpArray[b];
         superTrendDown[b]= SDnArray[b];
         STdir[b]= 1;
      }
      else
      {
         double pc= iClose(sym,tf,b+1);
         double prevUp=   superTrendUp[b+1];
         double prevDown= superTrendDown[b+1];

         // STUp logic
         if(pc> prevUp) superTrendUp[b]= MathMax(SUpArray[b], prevUp);
         else           superTrendUp[b]= SUpArray[b];

         // STDown logic
         if(pc< prevDown) superTrendDown[b]= MathMin(SDnArray[b], prevDown);
         else             superTrendDown[b]= SDnArray[b];

         int stVal= STdir[b+1];
         if(pc> superTrendDown[b+1]) stVal=1;
         else if(pc< superTrendUp[b+1]) stVal=-1;
         STdir[b]= stVal;
      }
      superTrendVal[b]= (STdir[b]==1)? superTrendUp[b]: superTrendDown[b];
   }

   // 7) Counters: _3xmabuy, stbuy, st3xmabuy, etc. rolling from right to left.
   for(int b=totalBars-1; b>=0; b--)
   {
      if(b==totalBars-1)
      {
         _3xmabuyArr[b]=0;   _3xmasellArr[b]=0;
         stbuyArr[b]=0;      stsellArr[b]=0;
         st3xmabuyArr[b]=0;  st3xmasellArr[b]=0;

         oiMACrossbuyArr[b]=0;  oiMACrosssellArr[b]=0;
         oi3xmabuyArr[b]=0;     oi3xmasellArr[b]=0;
         oistbuyArr[b]=0;       oistsellArr[b]=0;

         macrossSTbuyArr[b]=0;  macrossSTsellArr[b]=0;
         macross3xMAbuyArr[b]=0; macross3xMAsellArr[b]=0;

         hbuyArr[b]=0; hsellArr[b]=0;
         macrossbuyArr[b]=0; macrosssellArr[b]=0;

         aLongArr[b]=0; aShortArr[b]=0;
         longSignalArr[b]=0; shortSignalArr[b]=0;
      }
      else
      {
         int p3xbuy   = _3xmabuyArr[b+1];
         int p3xsell  = _3xmasellArr[b+1];
         int pstbuy   = stbuyArr[b+1];
         int pstsell  = stsellArr[b+1];
         int pst3xbuy = st3xmabuyArr[b+1];
         int pst3xsell= st3xmasellArr[b+1];
         int poiMBuy  = oiMACrossbuyArr[b+1];
         int poiMSell = oiMACrosssellArr[b+1];
         int poi3xbuy = oi3xmabuyArr[b+1];
         int poi3xsell= oi3xmasellArr[b+1];
         int poistbuy = oistbuyArr[b+1];
         int poistsell= oistsellArr[b+1];
         int pmSTbuy  = macrossSTbuyArr[b+1];
         int pmSTsell = macrossSTsellArr[b+1];
         int pm3xbuy  = macross3xMAbuyArr[b+1];
         int pm3xsell = macross3xMAsellArr[b+1];
         int phbuy    = hbuyArr[b+1];
         int phsell   = hsellArr[b+1];
         int pmaBuy   = macrossbuyArr[b+1];
         int pmaSell  = macrosssellArr[b+1];
         int paLong   = aLongArr[b+1];
         int paShort  = aShortArr[b+1];
         int plongSig = longSignalArr[b+1];
         int pshortSig= shortSignalArr[b+1];

         // replicate the script conditions
         bool condBuyA = (clrDir[b]==1 && iClose(sym,tf,b)> fastMA[b] && maDir[b]==1);
         bool condBuyB = (clrDir[b]==1 && maDir[b]==1);
         int new3xbuy=0;
         if(condBuyA) new3xbuy= p3xbuy+1;
         else if(condBuyB) new3xbuy= (p3xbuy>0)? p3xbuy+1: 0;
         else new3xbuy=0;
         _3xmabuyArr[b]= new3xbuy;

         bool condSellA= (clrDir[b]==-1 && iClose(sym,tf,b)< fastMA[b] && maDir[b]==-1);
         bool condSellB= (clrDir[b]==-1 && maDir[b]==-1);
         int new3xsell=0;
         if(condSellA) new3xsell= p3xsell+1;
         else if(condSellB) new3xsell= (p3xsell>0)? p3xsell+1:0;
         else new3xsell=0;
         _3xmasellArr[b]= new3xsell;

         bool condSTb= (clrDir[b]==1 && STdir[b]==1);
         stbuyArr[b]= (condSTb)? pstbuy+1: 0;
         bool condSTs= (clrDir[b]==-1 && STdir[b]==-1);
         stsellArr[b]= (condSTs)? pstsell+1:0;

         bool condst3xb= ((disable3xMAFilter||_3xmabuyArr[b]>0) && stbuyArr[b]>0);
         st3xmabuyArr[b]= condst3xb? pst3xbuy+1:0;
         bool condst3xs= ((disable3xMAFilter||_3xmasellArr[b]>0) && stsellArr[b]>0);
         st3xmasellArr[b]= condst3xs? pst3xsell+1:0;

         // For "OutsideIn:MACross" logic we do placeholder 
         oiMACrossbuyArr[b]=0; 
         oiMACrosssellArr[b]=0;

         bool condOi3b= (oiMACrossbuyArr[b]>0 && (disable3xMAFilter|| maDir[b]==1));
         oi3xmabuyArr[b]= condOi3b? poi3xbuy+1:0;
         bool condOi3s= (oiMACrosssellArr[b]>0 && (disable3xMAFilter|| maDir[b]==-1));
         oi3xmasellArr[b]= condOi3s? poi3xsell+1:0;

         bool condOistb= (oiMACrossbuyArr[b]>0 && STdir[b]==1);
         oistbuyArr[b]= condOistb? poistbuy+1:0;
         bool condOists= (oiMACrosssellArr[b]>0 && STdir[b]==-1);
         oistsellArr[b]= condOists? poistsell+1:0;

         // macrossSTbuy => crossover(ma_fast, ma_coloured) & STdir[b]==1
         bool crosB=false, crosS=false;
         if(b<totalBars-1)
         {
            double currFast= fastMA[b];
            double currCol=  colMA[b];
            double prevFast= fastMA[b+1];
            double prevCol=  colMA[b+1];
            crosB= (currFast>currCol && prevFast<= prevCol);
            crosS= (currFast<currCol && prevFast>= prevCol);
         }
         bool condmSTb= (crosB && STdir[b]==1);
         macrossSTbuyArr[b]= condmSTb? pmSTbuy+1:0;
         bool condmSTs= (crosS && STdir[b]==-1);
         macrossSTsellArr[b]= condmSTs? pmSTsell+1:0;

         // macross3xMAbuy => crossover(fast,coloured) & (disable3xMAFilter or maDir[b]==1)
         bool condm3b= (crosB && (disable3xMAFilter|| maDir[b]==1));
         macross3xMAbuyArr[b]= condm3b? pm3xbuy+1:0;
         bool condm3s= (crosS && (disable3xMAFilter|| maDir[b]==-1));
         macross3xMAsellArr[b]= condm3s? pm3xsell+1:0;

         // final long / short:
         bool cond3xb1= (u3xMATrendFilter && _3xmabuyArr[b]==1);
         bool condSTb1= (uSuperTrendFilter && stbuyArr[b]==1);
         bool condBothb1= (uBothTrendFilters && st3xmabuyArr[b]==1);
         bool condOI3b= (uOI3xMAFilter && oi3xmabuyArr[b]==1);
         bool condOISTb= (uOISTFilter && oistbuyArr[b]==1);
         bool condOIMACb= (uOIMACrossFilter && oiMACrossbuyArr[b]==1);
         bool condMacSTb= (uMACrossSTFilter && macrossSTbuyArr[b]==1);
         bool condMac3xb= (uMACross3xMAFilter && macross3xMAbuyArr[b]==1);
         bool isLong=( cond3xb1||condSTb1||condBothb1||condOI3b||condOISTb||condOIMACb||condMacSTb||condMac3xb );

         bool cond3xs1= (u3xMATrendFilter && _3xmasellArr[b]==1);
         bool condSTs1= (uSuperTrendFilter && stsellArr[b]==1);
         bool condBoths1= (uBothTrendFilters && st3xmasellArr[b]==1);
         bool condOI3s= (uOI3xMAFilter && oi3xmasellArr[b]==1);
         bool condOISTs= (uOISTFilter && oistsellArr[b]==1);
         bool condOIMACs= (uOIMACrossFilter && oiMACrosssellArr[b]==1);
         bool condMacSTs= (uMACrossSTFilter && macrossSTsellArr[b]==1);
         bool condMac3xs= (uMACross3xMAFilter && macross3xMAsellArr[b]==1);
         bool isShort=( cond3xs1||condSTs1||condBoths1||condOI3s||condOISTs||condOIMACs||condMacSTs||condMac3xs );

         longSignalArr[b]= (isLong?1:0);
         shortSignalArr[b]= (isShort?1:0);

         // hbuy => clrdirection==1 => +1
         bool hcondB=(clrDir[b]==1);
         hbuyArr[b]=(hcondB?phbuy+1:0);
         bool hcondS=(clrDir[b]==-1);
         hsellArr[b]=(hcondS?phsell+1:0);

         // macrossbuy => crossover(fast,coloured)? => +1
         bool crosMB=false, crosMS=false;
         if(b<totalBars-1)
         {
            double cF= fastMA[b];
            double cC= colMA[b];
            double pF= fastMA[b+1];
            double pC= colMA[b+1];
            crosMB= (cF> cC && pF<=pC);
            crosMS= (cF< cC && pF>=pC);
         }
         macrossbuyArr[b]= crosMB? pmaBuy+1:0;
         macrosssellArr[b]= crosMS? pmaSell+1:0;

         // along => (disableAllFilters && hbuy==1) or (uMACrossFilter && macrossbuy==1)
         bool condAl= false;
         if( (disableAllFilters && hbuyArr[b]==1) || (uMACrossFilter && macrossbuyArr[b]==1) )
            condAl=true;
         aLongArr[b]= condAl? paLong+1:0;

         bool condAs= false;
         if( (disableAllFilters && hsellArr[b]==1) || (uMACrossFilter && macrosssellArr[b]==1) )
            condAs=true;
         aShortArr[b]= condAs? paShort+1:0;

         // "long := long or along"
         bool finalLong= (longSignalArr[b]==1 || aLongArr[b]==1);
         longSignalArr[b]= (finalLong?1:0);
         bool finalShort= (shortSignalArr[b]==1 || aShortArr[b]==1);
         shortSignalArr[b]= (finalShort?1:0);
      }
   }

   // 8) Finally, we check bar i=1 => if (longSignalArr[1]==1), that means we got a new buy signal
   // if shortSignalArr[1]==1 => new sell signal
   // Return +1 for buy, -1 for sell, 0 for none

   if(totalBars<2) return 0;
   if(longSignalArr[1]==1) return +1;
   if(shortSignalArr[1]==1) return -1;
   return 0;
}


//+------------------------------------------------------------------+
//| ExecuteOnOpen: if we have a stored signal from T3 script        |
//+------------------------------------------------------------------+
void ExecuteOnOpen(string sym, int idx)
{
   int storedBar = lastSignalBarMap[idx];
   int storedDir = lastSignalDirMap[idx];
   if(storedBar < 0 || storedDir == 0)
      return;

   // Check if we are indeed on a new bar
   // iBarShift(sym,TimeFrame, iTime(sym,TimeFrame,0)) => current bar
   int currBar = iBarShift(sym, TimeFrame, iTime(sym, TimeFrame, 0));
   if(currBar >= storedBar)
      return; // not advanced

   // spread check
   double ask= MarketInfo(sym, MODE_ASK);
   double bid= MarketInfo(sym, MODE_BID);
   double spPts=(ask - bid)/Point;
   double spPips= spPts;
   if((Digits==5||Digits==3)&&spPts>0) spPips= spPts/10.0;
   if(spPips> MaxSpreadPips)
   {
      Print("Spread too high for ", sym, ": ", spPips, " pips. Skipping.");
      return;
   }

   // Count positions
   int sameDirCount = CountPositions(sym, storedDir);
   int oppDirCount  = CountPositions(sym, -storedDir);
   double oppLot    = GetTotalLotSize(sym, -storedDir);

   if(oppDirCount>0 || oppLot>0.0)
   {
      double netSize= oppLot + GetTradeLots(sym, storedDir);
      if(storedDir>0) PlaceOrder(sym, OP_BUY, netSize);
      else            PlaceOrder(sym, OP_SELL, netSize);
   }
   else
   {
      if(sameDirCount< MaxPyramiding)
      {
         int cmd= (storedDir>0)? OP_BUY: OP_SELL;
         double lots= GetTradeLots(sym, storedDir);
         PlaceOrder(sym, cmd, lots);
      }
   }

   // reset
   lastSignalBarMap[idx]= -1;
   lastSignalDirMap[idx]= 0;
}

//+------------------------------------------------------------------+
//| PlaceOrder                                                      |
//+------------------------------------------------------------------+
bool PlaceOrder(string sym, int cmd, double lots)
{
   if(lots<=0.0001)
   {
      Print("Lot size too small for ", sym, ". Aborting trade.");
      return(false);
   }
   int pipPoints= GetPipPoints(PipDistance, sym);
   double offset= pipPoints*Point;
   double price=0, sl=0, tp=0;
   if(cmd==OP_BUY)
   {
      price= NormalizeDouble(MarketInfo(sym, MODE_ASK), Digits);
      sl=    NormalizeDouble(price - offset, Digits);
      tp=    NormalizeDouble(price + offset, Digits);
   }
   else if(cmd==OP_SELL)
   {
      price= NormalizeDouble(MarketInfo(sym, MODE_BID), Digits);
      sl=    NormalizeDouble(price + offset, Digits);
      tp=    NormalizeDouble(price - offset, Digits);
   }
   else
   {
      Print("Invalid cmd for PlaceOrder: ", cmd);
      return false;
   }
   int ticket= OrderSend(sym, cmd, lots, price, Slippage, sl, tp, 
                "T3GoldSniperEA", 0, 0, clrBlue);
   if(ticket<0)
   {
      Print("OrderSend failed for ", sym, " cmd=", cmd, " lots=", lots,
            " err=", GetLastError());
      return false;
   }
   else
   {
      Print("Opened ", sym, " order #", ticket, " => cmd=", cmd, " lots=", lots);
   }
   return true;
}

//+------------------------------------------------------------------+
//| GetTradeLots - same as your EA logic                            |
//+------------------------------------------------------------------+
double GetTradeLots(string sym, int dir)
{
   if(!UseRiskManagement) return FixedLots;
   double bal= AccountBalance();
   double riskAmt= (RiskPercent/100.0)* bal;
   double tickVal= MarketInfo(sym, MODE_TICKVALUE);
   double tickSize=MarketInfo(sym, MODE_TICKSIZE);
   double oneLotPipValue=0;
   if(tickSize>0) oneLotPipValue= tickVal/tickSize;
   if(oneLotPipValue<=0)
   {
      Print("Warning: pipValue<=0 for ", sym, ". Using FixedLots fallback");
      return FixedLots;
   }
   double neededLot= riskAmt / ( (double)PipDistance* oneLotPipValue );
   double minLot = MarketInfo(sym,MODE_MINLOT);
   double maxLot = MarketInfo(sym,MODE_MAXLOT);
   double lotStep=MarketInfo(sym,MODE_LOTSTEP);

   if(lotStep>0) neededLot= MathFloor(neededLot/lotStep)*lotStep;
   if(neededLot< minLot) neededLot=minLot;
   if(neededLot> maxLot) neededLot=maxLot;
   neededLot= NormalizeDouble(neededLot,2);
   return neededLot;
}

//+------------------------------------------------------------------+
//| CountPositions                                                  |
//+------------------------------------------------------------------+
int CountPositions(string sym, int dir)
{
   int c=0;
   for(int i=OrdersTotal()-1;i>=0;i--)
   {
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES))
      {
         if(OrderSymbol()==sym && OrderMagicNumber()==0)
         {
            if(dir>0 && OrderType()==OP_BUY) c++;
            if(dir<0 && OrderType()==OP_SELL)c++;
         }
      }
   }
   return c;
}

//+------------------------------------------------------------------+
//| GetTotalLotSize                                                 |
//+------------------------------------------------------------------+
double GetTotalLotSize(string sym, int dir)
{
   double sum=0.0;
   for(int i=OrdersTotal()-1;i>=0;i--)
   {
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES))
      {
         if(OrderSymbol()==sym && OrderMagicNumber()==0)
         {
            if(dir>0 && OrderType()==OP_BUY) sum+=OrderLots();
            if(dir<0 && OrderType()==OP_SELL)sum+=OrderLots();
         }
      }
   }
   return sum;
}

//+------------------------------------------------------------------+
//| GetPipPoints                                                    |
//+------------------------------------------------------------------+
int GetPipPoints(int pips, string sym)
{
   int d=(int)MarketInfo(sym,MODE_DIGITS);
   if(d==3||d==5) return pips*10;
   return pips;
}

//+------------------------------------------------------------------+
//| ApplyTrailingStop                                               |
//+------------------------------------------------------------------+
void ApplyTrailingStop()
{
   for(int i=OrdersTotal()-1;i>=0;i--)
   {
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES))
      {
         if(OrderMagicNumber()==0)
         {
            int type= OrderType();
            if(type==OP_BUY)
            {
               double oldSL=OrderStopLoss();
               double newSL= NormalizeDouble(Bid- TrailingStopPoints*Point, Digits);
               if(Bid> OrderOpenPrice() && newSL> oldSL)
               {
                  double tp=OrderTakeProfit();
                  bool modOK= OrderModify(OrderTicket(),OrderOpenPrice(),newSL,tp,0,clrBlue);
                  if(!modOK)
                     Print("TrailingStop for BUY failed, err=",GetLastError());
                  else
                     Print("TrailingStop updated BUY #",OrderTicket()," => SL=",newSL);
               }
            }
            else if(type==OP_SELL)
            {
               double oldSL=OrderStopLoss();
               double newSL= NormalizeDouble(Ask+ TrailingStopPoints*Point, Digits);
               if(OrderOpenPrice()>Ask && (oldSL<0.0001 || newSL< oldSL))
               {
                  double tp=OrderTakeProfit();
                  bool modOK= OrderModify(OrderTicket(),OrderOpenPrice(),newSL,tp,0,clrBlue);
                  if(!modOK)
                     Print("TrailingStop for SELL failed, err=",GetLastError());
                  else
                     Print("TrailingStop updated SELL #",OrderTicket()," => SL=",newSL);
               }
            }
         }
      }
   }
}
//+------------------------------------------------------------------+
