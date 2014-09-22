/*
 * cloneable.h
 *
 *  Created on: Aug 18, 2014
 *      Author: ada
 */

#ifndef CLONEABLE_H_
#define CLONEABLE_H_

namespace SHEsis {
class Cloneable {
 public:
  virtual ~Cloneable() {}
  virtual Cloneable* clone(void) const = 0;
};

#define CLONE(X) \
  virtual X* clone(void) const { return new X(*this); }
}

#endif /* CLONEABLE_H_ */
