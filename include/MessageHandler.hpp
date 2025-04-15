/*
 *	This file is part of LCQPow.
 *
 *	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2022 by Jonas Hall et al.
 *
 *	LCQPow is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	LCQPow is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with LCQPow; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef LCQPOW_MESSAGEHANDLER_HPP
#define LCQPOW_MESSAGEHANDLER_HPP

#include "Utilities.hpp"

namespace LCQPow
{

class MessageHandler
{

public:
  /** Print a message.
   *
   * @param ret A return value which will impact the print message.
   * @param type The message type (will determine the message pre-print).
   *
   * @returns Simply passes the return value that was given. */
  static ReturnValue PrintMessage(ReturnValue ret, MessageType type = MESSAGE);

  /** Print a solution status.
   *
   * @param algoStat The status to be printed.
   *
   * @returns Simply passes the algorithm status that was given. */
  static AlgorithmStatus PrintSolution(AlgorithmStatus algoStat);
};
} // namespace LCQPow

#endif // LCQPOW_MESSAGEHANDLER_HPP
