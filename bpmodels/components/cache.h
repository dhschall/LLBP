/* MIT License
 *
 * Copyright (c) 2024 David Schall and EASE lab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */


template <typename key_t, typename value_t>
class BaseCache {
   protected:
    typedef typename std::pair<key_t, value_t> key_value_pair_t;
    typedef typename std::list<key_value_pair_t>::iterator list_iterator_t;
    typedef typename std::list<key_value_pair_t> set_t;

    std::unordered_map<uint64_t, list_iterator_t> _index;
    std::vector<std::list<key_value_pair_t>> _cache;
    const size_t _max_size;
    const size_t _assoc;
    const uint64_t _sets;
    const uint64_t _set_mask;

   public:
    BaseCache(size_t max_size, size_t assoc)
        : _max_size(max_size),
          _assoc(assoc),
          _sets(max_size / assoc),
          _set_mask(_sets - 1) {
        // Check if number of sets is a power of 2
        assert((_sets & (_sets - 1)) == 0);
        assert(_assoc * _sets == _max_size);
        _cache.resize(_sets);
        // for (auto& set : _cache) {
        //     set.resize(assoc);
        // }
    }

    void printCfg() {
        printf("Max size: %lu, Assoc: %lu, Sets: %lu\n", _max_size, _assoc,
               _sets);
    }

    size_t size() const { return _index.size(); }

    key_t index(const key_t& key) { return key & _set_mask; }

    set_t& getSet(const key_t& key) {
        return _cache[index(key)];
    }

    const std::unordered_map<key_t, list_iterator_t> getMap() { return _index; }

    value_t* get(const key_t& key) {
        auto it = _index.find(key);
        if (it == _index.end()) {
            return nullptr;
        }
        return &it->second->second;
    }

    void erase(const key_t& key) {
        auto it = _index.find(key);
        if (it == _index.end()) {
            return;
        }
        auto& set = getSet(key);
        set.erase(it->second);
        _index.erase(key);
    }

    value_t* getVictim(const key_t& key) {
        auto& set = getSet(key);
        if (set.size() < _assoc) {
            return nullptr;
        }
        return &set.back().second;
    }

    void touch(const key_t& key) {
        auto it = _index.find(key);
        if (it == _index.end()) {
            return;
        }
        auto& set = getSet(key);
        set.splice(set.begin(), set, it->second);
    }

    bool exists(const key_t& key) const {
        return _index.find(key) != _index.end();
    }

    int distance(const key_t& key) {
        auto it = _index.find(key);
        if (it == _index.end()) {
            return -1;
        }
        auto& set = getSet(key);
        return std::distance(set.begin(), it->second);
    }

    set_t& getResizedSet(const key_t& key) {
        auto& set = getSet(key);

        // If this element will exceed the max size, remove the last element
        if (set.size() >= _assoc) {
            auto last = set.end();
            last--;
            _index.erase(last->first);
            set.pop_back();
        }
        return set;
    }

    value_t* insertAt(const key_t& key, int at = 0) {
        auto v = get(key);
        if (v != nullptr) {
            return v;
        }

        // Get the set with a free item
        auto& set = getResizedSet(key);

        // Move to the insert position
        auto it2 = set.begin();
        at = std::min(at, (int)set.size());
        std::advance(it2, at);

        it2 = set.emplace(it2, key_value_pair_t(key, value_t()));
        _index[key] = it2;
        return &(it2->second);
    }

    value_t* insert(const key_t& key) {
        auto v = get(key);
        if (v != nullptr) {
            return v;
        }

        // Get the set with a free item
        auto& set = getResizedSet(key);

        // Move to the insert position
        auto it = set.begin();

        it = set.emplace(it, key_value_pair_t(key, value_t()));
        _index[key] = it;
        return &(it->second);
    }
};
