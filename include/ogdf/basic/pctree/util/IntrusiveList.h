#pragma once

#include <cassert>
#include <cstddef>

namespace pc_tree {
    template<class T>
    class IntrusiveList {
        T* first = nullptr;
        T* last = nullptr;
        size_t count = 0;

    public:
        class iterator {
            T* node;

        public:
            // iterator traits
            using iterator_category = std::input_iterator_tag;
            using value_type = T*;
            using difference_type = std::ptrdiff_t;
            using pointer = const T*;
            using reference = T*;

            explicit iterator(T* p_node) : node(p_node) {
            }

            iterator& operator++() {
                node = node->next;
                return *this;
            }

            iterator operator++(int) {
                iterator other = *this;
                ++(*this);
                return other;
            }

            bool operator==(iterator other) const {
                return node == other.node;
            }

            bool operator!=(iterator other) const {
                return node != other.node;
            }

            T* operator*() const {
                return node;
            }
        };

        class node {
            friend class IntrusiveList;

        private:
            T* next;
            T* prev;
        };

        iterator begin() const {
            return iterator(first);
        }

        iterator end() const {
            return iterator(nullptr);
        }

        void clear() {
            first = nullptr;
            last = nullptr;
            count = 0;
        }

        [[nodiscard]] bool empty() const {
            return count == 0;
        }

        [[nodiscard]] size_t size() const {
            return count;
        }

        T* front() const {
            assert(first != nullptr);
            return first;
        }

        T* back() const {
            assert(last != nullptr);
            return last;
        }

        void push_front(T* obj) {
            assert(obj != nullptr);
            check();

            if (first == nullptr) {
                obj->next = nullptr;
                last = obj;
            } else {
                obj->next = first;
                first->prev = obj;
            }

            obj->prev = nullptr;
            first = obj;
            count++;
            check();
        }

        void push_back(T* obj) {
            assert(obj != nullptr);
            check();

            if (last == nullptr) {
                obj->prev = nullptr;
                first = obj;
            } else {
                obj->prev = last;
                last->next = obj;
            }

            obj->next = nullptr;
            last = obj;
            count++;
            check();
        }

        void pop_front() {
            erase(front());
        }

        void pop_back() {
            erase(back());
        }

        void erase(T* obj) {
            assert(obj != nullptr);
            assert(count > 0);
            check();

            if (obj == first)
                first = obj->next;

            if (obj == last)
                last = obj->prev;

            if (obj->prev)
                obj->prev->next = obj->next;

            if (obj->next)
                obj->next->prev = obj->prev;

            count--;
            check();
        }

        void splice(iterator at, IntrusiveList<T>& other) {
            assert(at == begin() || at == end());
            check();
            other.check();

            if (at == end()) {
                if (last != nullptr) {
                    last->next = other.first;
                    other.first->prev = last;
                    last = other.last;
                } else {
                    first = other.first;
                    last = other.last;
                }
            } else if (at == begin()) {
                if (first != nullptr) {
                    first->prev = other.last;
                    other.last->next = first;
                    first = other.first;
                } else {
                    first = other.first;
                    last = other.last;
                }
            }

            count += other.count;
            other.count = 0;
            other.first = nullptr;
            other.last = nullptr;
            check();
            other.check();
        }

    private:
        void check() {
#ifndef NDEBUG
            size_t counter = 0;
            for ([[maybe_unused]] T* _ : *this)
                counter++;

            assert(counter == counter);
#endif
        }
    };
}