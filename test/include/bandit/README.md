bandit
======
[![Build Status](https://travis-ci.org/joakimkarlsson/bandit.png)](https://travis-ci.org/joakimkarlsson/bandit)

Human friendly unit testing for C++11

Bandit is a framework for C++11 that wants to make working with unit tests a pleasant
experience.

For more information, go to [the bandit website](http://banditcpp.org).

Bandit is released under the [MIT license](LICENSE.md)

#An example

This is a complete test application written in bandit:

```cpp
#include <bandit/bandit.h>
using namespace bandit;

// Tell bandit there are tests here.
go_bandit([](){

    // We're describing how a fuzzbox works.
    describe("fuzzbox:", [](){
      guitar_ptr guitar;
      fuzzbox_ptr fuzzbox;

      // Make sure each test has a fresh setup with
      // a guitar with a fuzzbox connected to it.
      before_each([&](){
        guitar = guitar_ptr(new struct guitar());
        fuzzbox = fuzzbox_ptr(new struct fuzzbox());
        guitar->add_effect(*fuzzbox);
      });

      it("starts in clean mode", [&](){
        AssertThat(guitar->sound(), Equals(sounds::clean));
      });

      // Describe what happens when we turn on the fuzzbox.
      describe("in distorted mode", [&](){

        // Turn on the fuzzbox.
        before_each([&](){
          fuzzbox->flip();
        });

        it("sounds distorted", [&](){
          AssertThat(guitar->sound(), Equals(sounds::distorted));
        });
      });
    });

});

int main(int argc, char* argv[])
{
  // Run the tests.
  return bandit::run(argc, argv);
}
```
