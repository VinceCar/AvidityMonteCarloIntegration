#include <iostream>
#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>

class ThreadPool {
public:
    // Constructor to initialize the thread pool with a fixed number of threads
    ThreadPool(size_t numThreads) {
        start(numThreads);
    }

    // Destructor to clean up threads
    ~ThreadPool() {
        stop();
    }

    // Submit a task to the pool
    template<class T>
    auto enqueue(T task) -> std::future<decltype(task())> {
        // Create a packaged task so we can retrieve the result via future
        auto wrapper = std::make_shared<std::packaged_task<decltype(task())()>>(std::move(task));

        {
            // Lock the queue and push the task
            std::unique_lock<std::mutex> lock(mEventMutex);
            mTasks.emplace([=] {
                (*wrapper)();
            });
        }

        // Notify one worker thread that a new task is available
        mEventVar.notify_one();
        
        // Return the future associated with the task
        return wrapper->get_future();
    }

private:
    // Vector of worker threads
    std::vector<std::thread> mThreads;

    // Task queue
    std::queue<std::function<void()>> mTasks;

    // Synchronization
    std::mutex mEventMutex;
    std::condition_variable mEventVar;

    // Flag to stop threads
    bool mStopping = false;

    // Start the pool with a fixed number of threads
    void start(size_t numThreads) {
        for (size_t i = 0; i < numThreads; ++i) {
            mThreads.emplace_back([=] {
                while (true) {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(mEventMutex);

                        // Wait for tasks to be available or for the pool to stop
                        mEventVar.wait(lock, [=] { return mStopping || !mTasks.empty(); });

                        if (mStopping && mTasks.empty()) {
                            break;
                        }

                        task = std::move(mTasks.front());
                        mTasks.pop();
                    }

                    // Execute the task
                    task();
                }
            });
        }
    }

    // Stop all threads and clean up
    void stop() {
        {
            std::unique_lock<std::mutex> lock(mEventMutex);
            mStopping = true;
        }

        mEventVar.notify_all();

        // Join all threads
        for (auto &thread : mThreads) {
            thread.join();
        }
    }
};
