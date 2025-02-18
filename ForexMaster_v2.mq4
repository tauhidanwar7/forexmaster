//+------------------------------------------------------------------+
//|  ForexMaster_v4_0_Improved_WithTrailingStop.mq4                 |
//|  MQL4 code with trailing stop (50 points).                      |
//+------------------------------------------------------------------+
#property strict

//--- Strategy Inputs
extern string  SymbolsToTrade        = "USDJPY,GBPJPY,EURJPY,GBPUSD";  
extern int     TimeFrame             = PERIOD_M30;

extern int     BB_Period             = 30;
extern double  BB_StDev              = 1.5;

extern int     ADX_Lookback          = 65;
extern int     ADX_EmaFast           = 4;
extern int     ADX_EmaSlow           = 8;

//--- Money Management
extern bool    UseRiskManagement     = false;        
extern double  RiskPercent           = 1.0;         
extern double  FixedLots             = 0.04;        
extern int     MaxPyramiding         = 5;           

//--- Trade Settings
extern int     PipDistance           = 60;
extern double  MaxSpreadPips         = 3.0;         
extern int     Slippage              = 3;           

//--- Trailing Stop Settings
extern bool    UseTrailingStop       = false;
extern int     TrailingStopPoints    = 50;

//--- Internal Tracking
datetime lastBarTimeMap[256];
int      lastSignalBarMap[256];
int      lastSignalDirMap[256];

//--- Symbol list array
string symbolList[];

//+------------------------------------------------------------------+
//| A replacement for SymbolInfoExists() in MQL4                    |
//| Checks if the symbol is recognized by the broker.               |
//+------------------------------------------------------------------+
bool MySymbolExists(string sym)
{
   ResetLastError();
   double digits = MarketInfo(sym, MODE_DIGITS);
   // error 4106 => unknown symbol
   if(GetLastError() == 4106)
      return false;
   return true;
}

//+------------------------------------------------------------------+
//| init() - EA initialization                                      |
//+------------------------------------------------------------------+
int init()
{
   // Parse the comma-separated symbol list
   StringSplit(SymbolsToTrade, ',', symbolList);

   // Initialize tracking arrays
   for(int i=0; i<ArraySize(symbolList); i++)
   {
      lastBarTimeMap[i]   = 0;
      lastSignalBarMap[i] = -1;
      lastSignalDirMap[i] =  0;
   }
   
   return(INIT_SUCCEEDED);
}

//+------------------------------------------------------------------+
//| deinit() - Cleanup                                              |
//+------------------------------------------------------------------+
int deinit()
{
   return(0);
}

//+------------------------------------------------------------------+
//| start() - Main function, called on every incoming tick          |
//+------------------------------------------------------------------+
int start()
{
   // For each symbol in our list, run the strategy logic
   for(int i=0; i<ArraySize(symbolList); i++)
   {
      string sym = symbolList[i];
      // Check if the symbol actually exists in MQL4
      if(MySymbolExists(sym))
      {
         ProcessSymbol(sym, i);
      }
   }

   // After processing all symbols, apply trailing stop if enabled
   if(UseTrailingStop)
      ApplyTrailingStop();

   return(0);
}

//+------------------------------------------------------------------+
//| ProcessSymbol: coordinate logic for each individual symbol      |
//+------------------------------------------------------------------+
void ProcessSymbol(string sym, int idx)
{
   // 1) Detect if there is a new bar on 'sym' (TimeFrame)
   datetime currentBarTime = iTime(sym, TimeFrame, 0);
   if(currentBarTime <= lastBarTimeMap[idx])
      return;  // no new bar or same bar
   
   // We have a new bar => update lastBarTime
   lastBarTimeMap[idx] = currentBarTime;

   // 2) Calculate the custom ADX buffers for the entire chart
   double adxFast[], adxSlow[];
   bool adxReady = CalculateCustomADX(sym, TimeFrame, adxFast, adxSlow);
   if(!adxReady) return;

   // 3) Check for signals on the just-closed bar (i=1)
   CheckForSignal(sym, idx, adxFast, adxSlow);

   // 4) If we have a stored signal, place the trade now at the new bar open
   ExecuteOnOpen(sym, idx);
}

//+------------------------------------------------------------------+
//| CalculateCustomADX: manually replicate the Pine Script approach |
//| Returns false if not enough bars                                |
//+------------------------------------------------------------------+
bool CalculateCustomADX(string sym, int tf, double& adxFast[], double& adxSlow[])
{
   int bars = iBars(sym, tf);
   if(bars <= ADX_Lookback+2) 
      return(false);

   // Prepare arrays
   ArrayResize(adxFast, bars);
   ArrayResize(adxSlow, bars);

   double SmoothedTR[], SmoothedDMPlus[], SmoothedDMMinus[];
   double DIPlus[], DIMinus[], DX[];
   ArrayResize(SmoothedTR,    bars);
   ArrayResize(SmoothedDMPlus,bars);
   ArrayResize(SmoothedDMMinus,bars);
   ArrayResize(DIPlus,        bars);
   ArrayResize(DIMinus,       bars);
   ArrayResize(DX,            bars);

   // For convenience, set them as series
   ArraySetAsSeries(adxFast,          true);
   ArraySetAsSeries(adxSlow,          true);
   ArraySetAsSeries(SmoothedTR,       true);
   ArraySetAsSeries(SmoothedDMPlus,   true);
   ArraySetAsSeries(SmoothedDMMinus,  true);
   ArraySetAsSeries(DIPlus,           true);
   ArraySetAsSeries(DIMinus,          true);
   ArraySetAsSeries(DX,               true);

   // Initialize
   for(int k=bars-1; k>=0; k--)
   {
      adxFast[k]=0; 
      adxSlow[k]=0;
      SmoothedTR[k]=0; 
      SmoothedDMPlus[k]=0; 
      SmoothedDMMinus[k]=0;
      DIPlus[k]=0; 
      DIMinus[k]=0; 
      DX[k]=0;
   }

   // Rightmost bar = oldest
   SmoothedTR[bars-1]      = 0;
   SmoothedDMPlus[bars-1]  = 0;
   SmoothedDMMinus[bars-1] = 0;
   DIPlus[bars-1]          = 0;
   DIMinus[bars-1]         = 0;
   DX[bars-1]              = 0;
   adxFast[bars-1]         = 0;
   adxSlow[bars-1]         = 0;

   // Loop from oldest to newest
   for(int i=bars-2; i>=0; i--)
   {
      double high_i   = iHigh(sym, tf, i);
      double low_i    = iLow(sym, tf, i);
      double close_i  = iClose(sym, tf, i);

      double high_i1  = iHigh(sym, tf, i+1);
      double low_i1   = iLow(sym, tf, i+1);
      double close_i1 = iClose(sym, tf, i+1);

      // True Range
      double tr1 = high_i - low_i;
      double tr2 = MathAbs(high_i - close_i1);
      double tr3 = MathAbs(low_i  - close_i1);
      double trueRange = MathMax(MathMax(tr1, tr2), tr3);

      // +DM / -DM
      double plusDM  = 0;
      double minusDM = 0;
      if((high_i - high_i1) > (low_i1 - low_i))
         plusDM = MathMax(high_i - high_i1, 0);
      if((low_i1 - low_i) > (high_i - high_i1))
         minusDM = MathMax(low_i1 - low_i, 0);

      double prevTR     = SmoothedTR[i+1];
      double prevDMPlus = SmoothedDMPlus[i+1];
      double prevDMMinus= SmoothedDMMinus[i+1];

      // Smooth them like in Pine Script:
      SmoothedTR[i]      = prevTR - (prevTR/ADX_Lookback) + trueRange;
      SmoothedDMPlus[i]  = prevDMPlus - (prevDMPlus/ADX_Lookback) + plusDM;
      SmoothedDMMinus[i] = prevDMMinus - (prevDMMinus/ADX_Lookback) + minusDM;

      // DI+ / DI-
      if(SmoothedTR[i] > 0)
      {
         DIPlus[i]  = (SmoothedDMPlus[i]  / SmoothedTR[i]) * 100.0;
         DIMinus[i] = (SmoothedDMMinus[i] / SmoothedTR[i]) * 100.0;
      }

      // DX
      double denom = DIPlus[i] + DIMinus[i];
      if(denom > 0)
         DX[i] = (MathAbs(DIPlus[i] - DIMinus[i]) / denom)*100.0;
      else
         DX[i] = 0;

      // ADX (Fast / Slow) => EMA of DX
      if(i == bars-2)  // initialize
      {
         adxFast[i] = DX[i];
         adxSlow[i] = DX[i];
      }
      else
      {
         adxFast[i] = EmaOnValue(adxFast[i+1], DX[i], ADX_EmaFast);
         adxSlow[i] = EmaOnValue(adxSlow[i+1], DX[i], ADX_EmaSlow);
      }
   }

   return(true);
}

//+------------------------------------------------------------------+
//| EmaOnValue: incremental EMA                                     |
//+------------------------------------------------------------------+
double EmaOnValue(double prevEma, double value, int period)
{
   double alpha = 2.0 / (period+1.0);
   return (prevEma + alpha*(value - prevEma));
}

//+------------------------------------------------------------------+
//| CheckForSignal: on bar i=1, see if cross above/below Bollinger  |
//+------------------------------------------------------------------+
void CheckForSignal(string sym, int idx, double& adxFast[], double& adxSlow[])
{
   // We use the just-closed bar: i=1
   int i = 1;
   int bars = iBars(sym, TimeFrame);
   if(i>=bars) return;

   // Price at bar i
   double close_i = iClose(sym, TimeFrame, i);

   // Bollinger logic
   double sma   = iMA(sym, TimeFrame, BB_Period, 0, MODE_SMA, PRICE_CLOSE, i);
   double stdev = iStdDev(sym, TimeFrame, BB_Period, 0, MODE_SMA, PRICE_CLOSE, i);
   double lowerBand = sma - (BB_StDev * stdev);
   double upperBand = sma + (BB_StDev * stdev);

   // ADX filter => only trade if ADX_Fast < ADX_Slow
   double f = adxFast[i];
   double s = adxSlow[i];
   bool adxOk = (f < s);

   // Check crossOver / crossUnder
   double prevClose = iClose(sym, TimeFrame, i+1);
   double prevSMA   = iMA(sym, TimeFrame, BB_Period, 0, MODE_SMA, PRICE_CLOSE, i+1);
   double prevStd   = iStdDev(sym, TimeFrame, BB_Period, 0, MODE_SMA, PRICE_CLOSE, i+1);
   double prevLB    = prevSMA - (BB_StDev * prevStd);
   double prevUB    = prevSMA + (BB_StDev * prevStd);

   bool crossOver  = (prevClose <= prevLB && close_i > lowerBand);
   bool crossUnder = (prevClose >= prevUB && close_i < upperBand);

   int newSignal = 0; // +1=Long, -1=Short
   if(crossOver && adxOk)
      newSignal = +1;
   else if(crossUnder && adxOk)
      newSignal = -1;

   if(newSignal != 0)
   {
      // Store signal for execution at next bar open
      lastSignalBarMap[idx] = iBarShift(sym, TimeFrame, iTime(sym, TimeFrame, i));
      lastSignalDirMap[idx] = newSignal;
   }
}

//+------------------------------------------------------------------+
//| ExecuteOnOpen: if we have a stored signal, trade at new bar open|
//+------------------------------------------------------------------+
void ExecuteOnOpen(string sym, int idx)
{
   int storedBar = lastSignalBarMap[idx];
   int storedDir = lastSignalDirMap[idx];
   if(storedBar < 0 || storedDir == 0)
      return;

   // Check if we are indeed on a new bar
   int currBar = iBarShift(sym, TimeFrame, iTime(sym, TimeFrame, 0));
   // If currBar >= storedBar, we haven't advanced
   if(currBar >= storedBar)
      return;

   // Spread check (skip trades if spread > MaxSpreadPips)
   double ask = MarketInfo(sym, MODE_ASK);
   double bid = MarketInfo(sym, MODE_BID);
   double spreadPoints = (ask - bid)/Point; 
   double spreadPips   = spreadPoints;
   // For 5-digit or 3-digit quotes, 10 points = 1 pip
   if((Digits == 5 || Digits == 3) && spreadPoints > 0) 
      spreadPips = spreadPoints / 10.0;
   if(spreadPips > MaxSpreadPips)
   {
      Print("Spread too high for ", sym, ": ", spreadPips, " pips. Skipping.");
      return;
   }

   // Opposite vs. same direction positions
   int sameDirCount = CountPositions(sym, storedDir);
   int oppDirCount  = CountPositions(sym, -storedDir);
   double sameLot   = GetTotalLotSize(sym, storedDir);
   double oppLot    = GetTotalLotSize(sym, -storedDir);

   // If there are positions in the opposite direction, net them out
   if(oppDirCount>0 || oppLot>0.0)
   {
      double netSize = oppLot + GetTradeLots(sym, storedDir);
      // If new signal is Long, we place a Buy netSize
      // If new signal is Short, we place a Sell netSize
      if(storedDir > 0)
         PlaceOrder(sym, OP_BUY, netSize);
      else
         PlaceOrder(sym, OP_SELL, netSize);
   }
   else
   {
      // If no opposite positions, check if we can pyramid
      if(sameDirCount < MaxPyramiding)
      {
         int cmd = (storedDir>0) ? OP_BUY : OP_SELL;
         double lots = GetTradeLots(sym, storedDir);
         PlaceOrder(sym, cmd, lots);
      }
   }

   // Reset the stored signal
   lastSignalBarMap[idx] = -1;
   lastSignalDirMap[idx] =  0;
}

//+------------------------------------------------------------------+
//| PlaceOrder: single order with PipDistance bracket SL & TP            |
//+------------------------------------------------------------------+
bool PlaceOrder(string sym, int cmd, double lots)
{
   if(lots <= 0.0001)
   {
      Print("Lot size too small for ", sym, ". Aborting trade.");
      return(false);
   }

   // Convert 50 pips to points
   int pipPoints = GetPipPoints(PipDistance, sym);
   double offset = pipPoints * Point;

   double price=0, sl=0, tp=0;
   if(cmd == OP_BUY)
   {
      price = NormalizeDouble(MarketInfo(sym, MODE_ASK), Digits);
      sl    = NormalizeDouble(price - offset, Digits);
      tp    = NormalizeDouble(price + offset, Digits);
   }
   else if(cmd == OP_SELL)
   {
      price = NormalizeDouble(MarketInfo(sym, MODE_BID), Digits);
      sl    = NormalizeDouble(price + offset, Digits);
      tp    = NormalizeDouble(price - offset, Digits);
   }
   else
   {
      Print("Invalid cmd for PlaceOrder: ", cmd);
      return false;
   }

   int ticket = OrderSend(sym, cmd, lots, price, Slippage, sl, tp,
                          "ForexMaster_v4.0_Improved", 0, 0, clrBlue);
   if(ticket<0)
   {
      Print("OrderSend failed for ", sym, " (cmd=", cmd, ", lots=", lots, 
            ") err=", GetLastError());
      return false;
   }
   else
   {
      Print("Opened ", sym, " order #", ticket, ": cmd=", cmd, " lots=", lots);
   }

   return true;
}

//+------------------------------------------------------------------+
//| GetTradeLots: either risk-based or fixed-lot                    |
//+------------------------------------------------------------------+
double GetTradeLots(string sym, int dir)
{
   if(!UseRiskManagement)
      return FixedLots;  // user-defined fixed lots

   // If we are using risk-based sizing:
   // - We assume a 50-pip stop for the bracket
   // - Risk = (RiskPercent% of balance)
   // => lotSize = (riskAmount) / (Value per pip * 50 pips)
   double balance = AccountBalance();
   double riskAmount = (RiskPercent/100.0) * balance;
   
   // For 1 lot, the pip value in account currency:
   double tickValue = MarketInfo(sym, MODE_TICKVALUE);
   double tickSize  = MarketInfo(sym, MODE_TICKSIZE);
   double oneLotPipValue = 0;
   if(tickSize > 0.0)
      oneLotPipValue = tickValue / tickSize;
   else
      oneLotPipValue = 0.0;

   if(oneLotPipValue <= 0)
   {
      // fallback if we can't calculate pip value
      Print("Warning: pipValue <= 0 for ", sym, ". Using FixedLots fallback.");
      return FixedLots;
   }

   // Our stop is 50 pips => risk for X lots = X * 50 * oneLotPipValue
   double neededLot = riskAmount / ( (double)PipDistance * oneLotPipValue );

   // Adjust for broker constraints
   double minLot  = MarketInfo(sym, MODE_MINLOT);
   double maxLot  = MarketInfo(sym, MODE_MAXLOT);
   double lotStep = MarketInfo(sym, MODE_LOTSTEP);

   // Round down to nearest lotStep
   if(lotStep > 0)
      neededLot = MathFloor(neededLot/lotStep) * lotStep;

   if(neededLot < minLot) neededLot = minLot;
   if(neededLot > maxLot) neededLot = maxLot;

   // Typically 2 decimal places for lot sizing
   neededLot = NormalizeDouble(neededLot, 2);

   return neededLot;
}

//+------------------------------------------------------------------+
//| CountPositions: number of open trades in given direction         |
//| dir=+1 => buys, dir=-1 => sells                                  |
//+------------------------------------------------------------------+
int CountPositions(string sym, int dir)
{
   int count=0;
   for(int i=OrdersTotal()-1; i>=0; i--)
   {
      if(OrderSelect(i, SELECT_BY_POS, MODE_TRADES))
      {
         if(OrderSymbol()==sym && OrderMagicNumber()==0)
         {
            if(dir>0 && OrderType()==OP_BUY)  count++;
            if(dir<0 && OrderType()==OP_SELL) count++;
         }
      }
   }
   return count;
}

//+------------------------------------------------------------------+
//| GetTotalLotSize: sum of open lots in a given direction           |
//+------------------------------------------------------------------+
double GetTotalLotSize(string sym, int dir)
{
   double total=0.0;
   for(int i=OrdersTotal()-1; i>=0; i--)
   {
      if(OrderSelect(i, SELECT_BY_POS, MODE_TRADES))
      {
         if(OrderSymbol()==sym && OrderMagicNumber()==0)
         {
            if(dir>0 && OrderType()==OP_BUY)  
               total += OrderLots();
            if(dir<0 && OrderType()==OP_SELL) 
               total += OrderLots();
         }
      }
   }
   return total;
}

//+------------------------------------------------------------------+
//| GetPipPoints: convert pipDistance to broker-specific points      |
//+------------------------------------------------------------------+
int GetPipPoints(int pips, string sym)
{
   // If symbol has 5 or 3 decimal digits, 1 pip = 10 points
   int d = (int)MarketInfo(sym, MODE_DIGITS);
   if(d == 3 || d == 5)
      return (pips*10);
   return pips;
}

//+------------------------------------------------------------------+
//| ApplyTrailingStop: loops through open trades and updates SL     |
//+------------------------------------------------------------------+
void ApplyTrailingStop()
{
   // For each open order, if in profit, tighten SL to trail price by TrailingStopPoints
   // This is a "simple" trailing stop: 
   //   For BUY: newStop = Bid - TrailingStopPoints*Point
   //   For SELL: newStop = Ask + TrailingStopPoints*Point
   // We only move the stop if it improves it (i.e., makes it bigger on a BUY, or smaller on a SELL).
   for(int i=OrdersTotal()-1; i>=0; i--)
   {
      if(OrderSelect(i, SELECT_BY_POS, MODE_TRADES))
      {
         if(OrderMagicNumber()==0)
         {
            int type   = OrderType();
            string sym = OrderSymbol();
            
            if(type == OP_BUY)
            {
               // current stop
               double oldSL = OrderStopLoss();
               // new suggested stop
               double newSL = NormalizeDouble(Bid - TrailingStopPoints*Point, Digits);

               // only move SL up if newSL > oldSL and the position is in profit
               // "position in profit" => Bid > OrderOpenPrice()
               if(Bid > OrderOpenPrice() && newSL > oldSL)
               {
                  // We keep TP the same
                  double tp = OrderTakeProfit();
                  bool modOK = OrderModify(OrderTicket(), OrderOpenPrice(), newSL, tp, 0, clrBlue);
                  if(!modOK)
                  {
                     Print("Trailing stop for BUY failed, err=", GetLastError());
                  }
                  else
                  {
                     Print("Trailing stop updated for BUY #", OrderTicket(), " => SL=", newSL);
                  }
               }
            }
            else if(type == OP_SELL)
            {
               double oldSL = OrderStopLoss();
               double newSL = NormalizeDouble(Ask + TrailingStopPoints*Point, Digits);

               // only move SL down if newSL < oldSL and the position is in profit
               // "position in profit" => OrderOpenPrice() > Ask
               if(OrderOpenPrice() > Ask && (oldSL < 0.0001 || newSL < oldSL))
               {
                  double tp = OrderTakeProfit();
                  bool modOK = OrderModify(OrderTicket(), OrderOpenPrice(), newSL, tp, 0, clrBlue);
                  if(!modOK)
                  {
                     Print("Trailing stop for SELL failed, err=", GetLastError());
                  }
                  else
                  {
                     Print("Trailing stop updated for SELL #", OrderTicket(), " => SL=", newSL);
                  }
               }
            }
         }
      }
   }
}
//+------------------------------------------------------------------+
