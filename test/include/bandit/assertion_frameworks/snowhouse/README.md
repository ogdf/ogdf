snowhouse
=========

An assertion library for C++

Snowhouse is a stand alone assertion framework for C++. It was originally 
developed as part of [Igloo](http://github.com/joakimkarlsson/igloo) and has
been extracted to be usable in other contexts.

## Usage

```C++
#include <snowhouse/snowhouse.h>
using namespace snowhouse;

int main()
{
  std::cout << "Testing that 23 is 23" << std::endl;
  AssertThat(23, Is().EqualTo(23));

  try
  {
    AssertThat(12, Is().LessThan(11).And().GreaterThan(99));
  }
  catch(const AssertionException& ex)
  {
    std::cout << "Apparently this failed:" << std::endl;
    std::cout << ex.GetMessage() << std::endl;
  }

  return 0;
}
```

### Assertions

Snowhouse uses a constraint based assertion model that is heavily inspired by the
model used in [NUnit](http://nunit.org/). An assertion in Snowhouse is written 
using the following format:

```cpp
AssertThat(actual_value, <constraint expression>);
```

where &lt;constraint expression&gt; is an expression that actual_value is evaluated against when the test is executed.

Constraint expressions come in two basic forms: composite and fluent expressions

#### Composite Expressions

With composite expressions, you can create compact, powerful expressions that combine a set of predefined constraints with ones that you provide yourself.

Example:

```cpp
AssertThat(length, IsGreaterThan(4) && !Equals(10));
```

Composite expressions can be any combination of constraints and the standard logical C++ operators.

You can also add your own constraints to be used within composite expressions.

####Fluent Expressions

With fluent expressions, you can create assertions that better convey the intent of a test without exposing implementation-specific details. Fluent expressions aim to help you create tests that are not just by developers for developers, but rather can be read and understood by domain experts.

Fluent expressions also has the ability to make assertions on the elements in a conteiner in a way you cannot achieve with composite expressions.

Example:

```cpp
AssertThat(length, Is().GreaterThan(4).And().Not().EqualTo(10));
```

### Basic Constraints

####Equality Constraint

Used to verify equality between actual and expected.

```cpp
AssertThat(x, Equals(12));
AssertThat(x, Is().EqualTo(12));
```

####EqualityWithDelta Constraint

Used to verify equality between actual and expected, allowing the two to differ by a delta.

```cpp
AssertThat(2.49, EqualsWithDelta(2.5, 0.1));
AssertThat(2.49, Is().EqualToWithDelta(2.5, 0.1));
```

####GreaterThan Constraint

Used to verify that actual is greater than a value.

```cpp
AssertThat(x, IsGreaterThan(4));
AssertThat(x, Is().GreaterThan(4));
```


####LessThan Constraint

Used to verify that actual is less than a value.

```cpp
AssertThat(x, IsLessThan(3));
AssertThat(x, Is().LessThan(3));
```

### String Constraints

String assertions in Snowhouse are used to verify the values of STL strings (std::string).

####Equality Constraints

Used to verify that actual is equal to an expected value.

```cpp
Assert:That(actual_str, Equals("foo")); 
AssertThat(actual_str, Is().EqualTo("foo"));
```

####Contains Constraint

Used to verify that a string contains a substring.

```cpp
AssertThat(actual_str, Contains("foo")); 
AssertThat(actual_str, Is().Containing("foo"));
```

####EndsWith Constraint

Used to verify that a string ends with an expected substring.

```cpp
AssertThat(actual_str, EndsWith("foo"));
AssertThat(actual_str, Is().EndingWith("foo"));
```

####StartsWith Constraint

Used to verify that a string starts with an expected substring.

```cpp
AssertThat(actual_str, StartsWith("foo"));
AssertThat(actual_str, Is().StartingWith("foo"));
```

####HasLength Constraint

Used to verify that a string is of the expected length.

```cpp
AssertThat(actual_str, HasLength(5));
AssertThat(actual_str, Is().OfLength(5));
```

###Constraints on Multi Line Strings

If you have a string that contains multiple lines, you can use the collection constraints to make assertions on the content of that string. This may be handy if you have a string that, for instance, represents the resulting content of a file or a network transmission.

Snowhouse can handle both windows (CR+LF) and unix (LF) line endings

```cpp
std::string lines = "First line\r\nSecond line\r\nThird line"; 
AssertThat(lines, Has().Exactly(1).StartingWith("Second"));
```

###Container Constraints

The following constraints can be applied to containers in the standard template library:

####Contains Constraint

Used to verify that a container contains an expected value.

```cpp
AssertThat(container, Contains(12));
AssertThat(container, Is().Containing(12));
```

####HasLength Constraint

Used to verify that a container has the expected length.

```cpp
AssertThat(container, HasLength(3));
AssertThat(container, Is().OfLength(3));
```

####IsEmpty Constraint

Used to verify that a container is empty.

```cpp
AssertThat(contatiner, IsEmpty());
AssertThat(container, Is().Empty());
```

####All

Used to verify that all elements of a STL sequence container matches an expectation.

```cpp
AssertThat(container, Has().All().LessThan(5).Or().EqualTo(66));
```

####AtLeast

Used to verify that at least a specified amount of elements in a STL sequence container matches an expectation.

```cpp
AssertThat(container, Has().AtLeast(3).StartingWith("foo"));
```

####AtMost

Used to verify that at most a specified amount of elements in a STL sequence container matches an expectation.

```cpp
Assert:That(container, Has().AtMost(2).Not().Containing("failed"));
```

####Exactly

Used to verify that a STL sequence container has exactly a specified amount of elements that matches an expectation.

```cpp
AssertThat(container, Has().Exactly(3).GreaterThan(10).And().LessThan(20));
```

####EqualsContainer

Used to verify that two STL sequence containers are equal.

```cpp
AssertThat(container1, EqualsContainer(container2));
AssertThat(container1, Is().EqualToContainer(container2));
```

#####Predicate functions

You can supply a predicate function or a functor to EqualsContainer to customize how to compare the elements in the two containers.

With a predicate function:

```cpp
static bool are_my_types_equal(const my_type& lhs, const my_type& rhs)
{
  return lhs.my_val_ == rhs.my_val_;
}

AssertThat(container1, EqualsContainer(container2, are_my_types_equal));
```

With a functor as predicate:

```cpp
struct within_delta
{
  within_delta(int delta) : delta_(delta) {}

  bool operator()(const my_type& lhs, const my_type& rhs) const
  {
    return abs(lhs.my_val_ - rhs.my_val_) <= delta_;
  }

private:
    int delta_;
};

AssertThat(container1, Is().EqualToContainer(container1, within_delta(1));
```

###Exceptions

Exception constraints can be used to verify that your code throws the correct exceptions.

####AssertThrows

AssertThrows succeeds if the exception thrown by the call is of the supplied type (or one of its subtypes).

```cpp
AssertThrows(std::logic_error, myObject.a_method(42));
```

####Making Assertions on the Thrown Exceptions

If AssertThrows succeeds, it will store the thrown exception so that you can make more detailed assertions on it.

```cpp
AssertThrows(std::logic_error, myObject.a_method(42));
AssertThat(LastException<std::logic_error>().what(), Is().Containing("logic failure"));
```

The LastException<> is available in the scope of the call to AssertThrows. An exception is not available between specs in order to avoid the result of one spec contaminating another.

###Custom Constraints

You can add your own constraints to Snowhouse to create more expressive specifications.

####Fulfills Constraints

By defining the following matcher

```cpp
struct IsEvenNumber
{
  bool Matches(const int actual) const
  {
    return (actual % 2) == 0; 
  }

  friend std::ostream& operator<<(std::ostream& stm, const IsEvenNumber& );
};

std::ostream& operator<<(std::ostream& stm, const IsEvenNumber& )
{
  stm << "An even number";
  return stm;
}
```

You can create the following constraints in Snowhouse:

```cpp
AssertThat(42, Fulfills(IsEvenNumber()));
AssertThat(42, Is().Fulfilling(IsEvenNumber()));
```

Your custom matcher should implement a method called Matches() that takes a parameter of the type you expect and returns true if the passed parameter fulfills the constraint.

To get more expressive failure messages, you should also implement the streaming operator as in the example above.
